"""Pool four recorder sites across the same four geometries as Figure 2B/2C."""

import itertools
import numpy as np
import pandas as pd

# Run from the repository root; reuse existing site mappings and sample annotations.
ROOT = "analysis/boxb_in_vitro_sequencing"
SELECTION_FILE = f"{ROOT}/annotations/fig2e_selection.csv"
SITE_FILE = f"{ROOT}/annotations/recorder_context_selection.csv"
POSITION_FILE = f"{ROOT}/annotations/recorder_positions.csv"
SAMPLE_FILE = f"{ROOT}/annotations/sample_info.csv"
BARCODE_FILE = f"{ROOT}/annotations/barcode_annotations.csv"
COUNT_DIRECTORY = f"{ROOT}/data/summary_stats_combined"
SITE_OUTPUT = f"{ROOT}/tables/fig2e_site_summary.csv"
OUTPUT_FILE = f"{ROOT}/tables/fig2e_context_summary.csv"
KEYS = ["sample_id", "position_id", "target_pos_to_boxb", "target_dist"]
EXPECTED_SAMPLES, EXPECTED_CONTEXTS, EXPECTED_BACKGROUNDS = 2, 16, 64
EXPECTED_SITES, INSERT_LENGTH = 4, 5
EXPECTED_GEOMETRIES = 4
REFERENCE_LIBRARY_TYPE = "boxb"
DNA_BASES, RNA_BASES = "ACGT", "ACGU"
READ_TO_RNA = dict(zip(DNA_BASES, "UGCA"))

selection = pd.read_csv(SELECTION_FILE)
if (len(selection) != EXPECTED_SAMPLES * EXPECTED_SITES * EXPECTED_GEOMETRIES
        or selection.sample_id.nunique() != EXPECTED_SAMPLES
        or selection.duplicated(KEYS).any()
        or not selection.groupby("sample_id").position_id.nunique().eq(EXPECTED_SITES).all()
        or len(selection[["target_pos_to_boxb", "target_dist"]].drop_duplicates()) != EXPECTED_GEOMETRIES
        or not selection.groupby(["sample_id", "position_id"]).size().eq(EXPECTED_GEOMETRIES).all()):
    raise ValueError("Select four sites at the same four geometries for two samples")
sites = pd.read_csv(SITE_FILE, dtype={"variable_subpos": str}).drop(columns=["target_pos_to_boxb", "target_dist"])
positions = pd.read_csv(POSITION_FILE).rename(columns={"context": "reference_context"})
choices = selection.merge(sites, on=["sample_id", "position_id"], how="left", validate="many_to_one")
choices = choices.merge(positions, on="position_id", how="left", validate="many_to_one")
if choices.isna().any().any():
    raise ValueError("Every selected site must resolve in the established context annotations")
samples = pd.read_csv(SAMPLE_FILE)
barcodes = pd.read_csv(BARCODE_FILE, dtype={"variable_subpos": str})
boxb_geometries = barcodes.loc[barcodes.variable_type.eq(REFERENCE_LIBRARY_TYPE),
                              ["target_pos_to_boxb", "target_dist", "target_length"]].drop_duplicates()
selected_geometries = selection[["target_pos_to_boxb", "target_dist"]].drop_duplicates()
if (len(boxb_geometries) != EXPECTED_GEOMETRIES
        or len(selected_geometries.merge(boxb_geometries, on=["target_pos_to_boxb", "target_dist"])) != EXPECTED_GEOMETRIES
        or boxb_geometries.target_length.nunique() != 1):
    raise ValueError("Selected geometries do not match the boxB-variation library")
frames = []
for sample_id, sample_choices in choices.groupby("sample_id", sort=False):
    sample = samples.loc[samples.sample_id.eq(sample_id)]
    if len(sample) != 1:
        raise ValueError("Sample metadata must resolve uniquely")
    for column in ["tada_type", "tada_conc", "condition"]:
        if not sample_choices[column].eq(sample.iloc[0][column]).all():
            raise ValueError(f"Sample metadata disagree: {column}")
    columns = ["sample_id", "barcode", "insert", "umi_counts"] + sorted(sample_choices.editing_column.unique())
    counts = pd.read_csv(f"{COUNT_DIRECTORY}/{sample_id}.csv.gz", usecols=columns)
    if not counts.sample_id.eq(sample_id).all():
        raise ValueError("Unexpected sample in count file")
    for choice in sample_choices.itertuples():
        barcode = barcodes.loc[barcodes.variable_type.eq(choice.variable_type)
                              & barcodes.variable_subpos.eq(choice.variable_subpos)
                              & barcodes.target_pos_to_boxb.eq(choice.target_pos_to_boxb)
                              & barcodes.target_dist.eq(choice.target_dist)]
        if len(barcode) != 1:
            raise ValueError("Each site and geometry must resolve to exactly one barcode")
        if barcode.iloc[0].target_length != boxb_geometries.target_length.iloc[0]:
            raise ValueError("Recorder length differs from the boxB-variation library")
        data = counts.loc[counts.barcode.eq(barcode.iloc[0].reverse_complement)].copy()
        if (data.empty or data["insert"].duplicated().any()
                or not data["insert"].str.fullmatch(f"[{DNA_BASES}]{{{INSERT_LENGTH}}}").all()
                or not data.umi_counts.gt(0).all()
                or not data[choice.editing_column].between(0, data.umi_counts).all()):
            raise ValueError("Invalid variant sequences or UMI counts")
        flank_indices = [choice.insert_fiveprime_index - 1, choice.insert_threeprime_index - 1]
        if len(set(flank_indices)) != 2 or not all(0 <= i < INSERT_LENGTH for i in flank_indices):
            raise ValueError("Invalid flanking-base indices")
        other_indices = [i for i in range(INSERT_LENGTH) if i not in flank_indices]
        backgrounds = {"".join(bases) for bases in itertools.product(DNA_BASES, repeat=len(other_indices))}
        data["rna_context"] = (data["insert"].str[flank_indices[0]].map(READ_TO_RNA)
                               + choice.reference_context[1]
                               + data["insert"].str[flank_indices[1]].map(READ_TO_RNA))
        data["background"] = data["insert"].map(lambda s: "".join(s[i] for i in other_indices))
        for context, rows in data.groupby("rna_context"):
            if len(rows) != EXPECTED_BACKGROUNDS or set(rows.background) != backgrounds:
                raise ValueError(f"{sample_id}, {choice.position_id}, {context} lacks 64 observed backgrounds")
        pooled = data.groupby("rna_context").agg(
            n_variants=("insert", "size"), total_umis=("umi_counts", "sum"),
            total_edited_umis=(choice.editing_column, "sum"),
            minimum_umis=("umi_counts", "min"), maximum_umis=("umi_counts", "max")).reset_index()
        for column in choices.columns:
            pooled[column] = getattr(choice, column)
        pooled["barcode"] = barcode.iloc[0].reverse_complement
        pooled["oligo_name"] = barcode.iloc[0].oligo_name
        pooled["pooled_fraction_central_a_edited"] = pooled.total_edited_umis / pooled.total_umis
        frames.append(pooled)
summary = pd.concat(frames, ignore_index=True)
if (summary.isna().any().any() or len(summary) != EXPECTED_SAMPLES * EXPECTED_SITES * EXPECTED_CONTEXTS * EXPECTED_GEOMETRIES
        or summary.duplicated(KEYS + ["rna_context"]).any()
        or not summary.n_variants.eq(EXPECTED_BACKGROUNDS).all()
        or len(summary[["tada_conc", "condition"]].drop_duplicates()) != 1):
    raise ValueError("Each sample must supply 16 contexts with 64 backgrounds at the same incubation")
if (not summary.total_umis.gt(0).all()
        or not summary.total_edited_umis.between(0, summary.total_umis).all()
        or not np.allclose(summary.pooled_fraction_central_a_edited,
                           summary.total_edited_umis / summary.total_umis, rtol=0, atol=1e-12)):
    raise ValueError("Editing fractions must match the pooled counts")
summary["fiveprime"] = summary.rna_context.str[0]
summary["threeprime"] = summary.rna_context.str[2]
for sample_id, data in summary.groupby("sample_id"):
    central = data.reference_context.str[1].unique()
    if len(central) != 1 or set(data.rna_context) != {
            f"{five}{central[0]}{three}" for five in RNA_BASES for three in RNA_BASES}:
        raise ValueError(f"{sample_id} does not contain all 16 contexts around the annotated central base")
summary.to_csv(SITE_OUTPUT, index=False)
summary["geometry"] = summary.target_pos_to_boxb.astype(str) + "prime_" + summary.target_dist.astype(str) + "nt"
metadata = ["sample_id", "tada_type", "tada_conc", "condition", "enzyme_label", "panel_order",
            "rna_context", "fiveprime", "threeprime"]
summary = summary.groupby(metadata, sort=False).agg(
    n_sites=("position_id", "nunique"),
    n_geometries=("geometry", "nunique"),
    geometries=("geometry", lambda values: ";".join(values.unique())),
    position_id=("position_id", lambda values: ";".join(values.unique())),
    reference_context=("reference_context", lambda values: ";".join(values.unique())),
    n_site_variants=("n_variants", "sum"), total_umis=("total_umis", "sum"),
    total_edited_umis=("total_edited_umis", "sum"),
    minimum_umis=("minimum_umis", "min"), maximum_umis=("maximum_umis", "max")
).reset_index()
if (len(summary) != EXPECTED_SAMPLES * EXPECTED_CONTEXTS or not summary.n_sites.eq(EXPECTED_SITES).all()
        or not summary.n_geometries.eq(EXPECTED_GEOMETRIES).all()
        or not summary.n_site_variants.eq(EXPECTED_SITES * EXPECTED_GEOMETRIES * EXPECTED_BACKGROUNDS).all()):
    raise ValueError("Every heat-map tile must include all four sites at all four geometries")
summary["pooled_fraction_central_a_edited"] = summary.total_edited_umis / summary.total_umis
summary = summary.sort_values(["panel_order", "rna_context"])
summary.to_csv(OUTPUT_FILE, index=False)
print(f"Pooled {EXPECTED_SITES} sites and {EXPECTED_GEOMETRIES} geometries into each of {len(summary)} context tiles")
print(summary.groupby("enzyme_label").pooled_fraction_central_a_edited.agg(["min", "max"]).to_string())
