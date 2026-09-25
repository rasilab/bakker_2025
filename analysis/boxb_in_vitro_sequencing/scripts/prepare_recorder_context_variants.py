"""Pool UMI counts across sequence backgrounds separately for each recorder site."""

import itertools
import pandas as pd

# Run from the repository root; experimental choices live in the selection CSV.
ROOT = "analysis/boxb_in_vitro_sequencing"
SELECTION_FILE = f"{ROOT}/annotations/recorder_context_selection.csv"
POSITION_FILE = f"{ROOT}/annotations/recorder_positions.csv"
SAMPLE_FILE = f"{ROOT}/annotations/sample_info.csv"
BARCODE_FILE = f"{ROOT}/annotations/barcode_annotations.csv"
COUNT_DIRECTORY = f"{ROOT}/data/summary_stats_combined"
VARIANT_FILE = f"{ROOT}/tables/recorder_context_variants.csv"
SUMMARY_FILE = f"{ROOT}/tables/recorder_context_summary.csv"
DNA_BASES, RNA_BASES = "ACGT", "ACGU"
READ_TO_RNA = dict(zip(DNA_BASES, "UGCA"))
INSERT_LENGTH, MINIMUM_UMIS, COVERAGE_REFERENCE = 5, 1, 50
EXPECTED_SAMPLES, EXPECTED_POSITIONS, EXPECTED_VARIANTS = 2, 4, 64

selection = pd.read_csv(SELECTION_FILE, dtype={"variable_subpos": str})
if (selection.sample_id.nunique() != EXPECTED_SAMPLES
        or len(selection) != EXPECTED_SAMPLES * EXPECTED_POSITIONS
        or selection.duplicated(["sample_id", "position_id"]).any()):
    raise ValueError("Select four distinct recorder sites for each of two samples")
shared_columns = ["tada_conc", "condition", "target_pos_to_boxb", "target_dist", "variable_type"]
if len(selection[shared_columns].drop_duplicates()) != 1:
    raise ValueError("Both enzymes must use the same concentration, incubation, and geometry")
site_columns = ["position_id", "variable_subpos", "insert_fiveprime_index", "insert_threeprime_index",
                "position_order", "marker"]
if len(selection[site_columns].drop_duplicates()) != EXPECTED_POSITIONS:
    raise ValueError("Both enzymes must use identical site definitions")
positions = pd.read_csv(POSITION_FILE).rename(columns={"context": "reference_context"})
selection = selection.merge(positions, on="position_id", how="left", validate="many_to_one")
if selection.editing_column.isna().any():
    raise ValueError("Every site must resolve in existing recorder annotations")
central_bases = selection.reference_context.str[1].unique()
if len(central_bases) != 1:
    raise ValueError("Selected sites do not share the same central base")
central_base = central_bases[0]
if (selection.reference_context.str[0].eq(central_base)
        | selection.reference_context.str[2].eq(central_base)).any():
    raise ValueError("A selected reference site has a neighboring central-base identity")
context_order = [f"{five}{central_base}{three}" for five, three in itertools.product(RNA_BASES, repeat=2)]
samples = pd.read_csv(SAMPLE_FILE)
annotations = pd.read_csv(BARCODE_FILE, dtype={"variable_subpos": str})
read_columns = ["sample_id", "barcode", "insert", "umi_counts"] + sorted(selection.editing_column.unique())
frames = []
for sample_id, choices in selection.groupby("sample_id", sort=False):
    sample = samples.loc[samples.sample_id.eq(sample_id)]
    if len(sample) != 1:
        raise ValueError("Sample must resolve uniquely in existing annotations")
    for column in ["tada_type", "tada_conc", "condition"]:
        if not choices[column].eq(sample.iloc[0][column]).all():
            raise ValueError(f"Sample metadata disagree with selection: {column}")
    counts = pd.read_csv(f"{COUNT_DIRECTORY}/{sample_id}.csv.gz", usecols=read_columns)
    if not counts.sample_id.eq(sample_id).all():
        raise ValueError("Count file contains an unexpected sample")
    for choice in choices.sort_values("position_order").itertuples():
        barcode_selection = annotations.loc[
            annotations.variable_type.eq(choice.variable_type)
            & annotations.variable_subpos.eq(choice.variable_subpos)
            & annotations.target_pos_to_boxb.eq(choice.target_pos_to_boxb)
            & annotations.target_dist.eq(choice.target_dist)]
        if len(barcode_selection) != 1:
            raise ValueError("Each site and geometry must resolve to exactly one barcode")
        barcode = barcode_selection.iloc[0]
        variants = counts.loc[counts.barcode.eq(barcode.reverse_complement),
                              ["sample_id", "barcode", "insert", "umi_counts", choice.editing_column]].copy()
        if variants.empty or variants["insert"].duplicated().any():
            raise ValueError("Expected one observed count row per unique insert")
        if not variants["insert"].str.fullmatch(f"[{DNA_BASES}]{{{INSERT_LENGTH}}}").all():
            raise ValueError("Noncanonical inserts require review before proceeding")
        if not variants.umi_counts.ge(MINIMUM_UMIS).all():
            raise ValueError("A variant lacks positive UMI coverage; review before proceeding")
        if not variants[choice.editing_column].between(0, variants.umi_counts).all():
            raise ValueError("Edited counts must be between zero and total UMI counts")
        flank_indices = [choice.insert_fiveprime_index - 1, choice.insert_threeprime_index - 1]
        if len(set(flank_indices)) != 2 or not all(0 <= index < INSERT_LENGTH for index in flank_indices):
            raise ValueError("Flanks must use two distinct insert positions")
        other_indices = [index for index in range(INSERT_LENGTH) if index not in flank_indices]
        expected_backgrounds = {"".join(bases) for bases in itertools.product(DNA_BASES, repeat=len(other_indices))}
        variants["fiveprime"] = variants["insert"].str[flank_indices[0]].map(READ_TO_RNA)
        variants["threeprime"] = variants["insert"].str[flank_indices[1]].map(READ_TO_RNA)
        variants["rna_context"] = variants.fiveprime + central_base + variants.threeprime
        variants["background_insert"] = variants["insert"].map(
            lambda insert: "".join(insert[index] for index in other_indices))
        if set(variants.rna_context) != set(context_order):
            raise ValueError("Not all 16 RNA contexts are present")
        for context, group in variants.groupby("rna_context"):
            if set(group.background_insert) != expected_backgrounds or len(group) != EXPECTED_VARIANTS:
                raise ValueError(f"{sample_id}, {choice.position_id}, {context} lacks 64 observed backgrounds")
        variants = variants.rename(columns={choice.editing_column: "edited_umi_counts"})
        variants["fraction_central_a_edited"] = variants.edited_umi_counts / variants.umi_counts
        variants["at_or_below_legacy_umi_cutoff"] = variants.umi_counts.le(COVERAGE_REFERENCE)
        for column in selection.columns:
            if column != "sample_id":
                variants[column] = getattr(choice, column)
        variants["oligo_name"] = barcode.oligo_name
        frames.append(variants)
        print(f"{choice.enzyme_label}, {choice.position_id}: {len(variants)} variants; "
              f"UMIs {variants.umi_counts.min()}–{variants.umi_counts.max()}")

variants = pd.concat(frames, ignore_index=True).sort_values(
    ["panel_order", "position_order", "rna_context", "background_insert"])
group_keys = ["sample_id", "position_id", "rna_context"]
summary = variants.groupby(group_keys, sort=True).agg(
    n_variants=("insert", "size"),
    minimum_umis=("umi_counts", "min"), median_umis=("umi_counts", "median"),
    maximum_umis=("umi_counts", "max"), total_umis=("umi_counts", "sum"),
    total_edited_umis=("edited_umi_counts", "sum"),
    n_at_or_below_legacy_umi_cutoff=("at_or_below_legacy_umi_cutoff", "sum"),
).reset_index()
summary["pooled_fraction_central_a_edited"] = summary.total_edited_umis / summary.total_umis
metadata_columns = list(selection.columns) + ["oligo_name"]
summary = summary.merge(variants[metadata_columns].drop_duplicates(),
                        on=["sample_id", "position_id"], validate="many_to_one")
summary = summary.sort_values(["panel_order", "position_order", "rna_context"])
variants.to_csv(VARIANT_FILE, index=False)
summary.to_csv(SUMMARY_FILE, index=False)
print(f"Wrote {len(summary)} pooled site/context measurements from {len(variants)} site/variant rows")
