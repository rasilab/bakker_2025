"""Prepare Recorder edit summary recorder and its time-course supplement; run from the repository root."""

import os

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control, ttest_rel

# Inputs and analysis choices. Preserve the published working baseline, including
# its 2–7-edit definition of the category labelled "2+ edits".
analysis_directory = "analysis/boxb_in_vitro_sequencing"
input_directory = f"{analysis_directory}/data/summary_stats_combined"
annotation_file = f"{analysis_directory}/annotations/barcode_annotations.csv"
sample_file = f"{analysis_directory}/annotations/sample_info.csv"
stem_file = f"{analysis_directory}/tables/boxb_wt_mut_stems.csv"
output_directory = f"{analysis_directory}/tables"
enzyme_order = ["lambdaN", "tada_only"]
concentration_order = ["125nM", "250nM", "500nM"]
insert_order = ["mut", "wt"]
edit_order = ["frac_1edit", "frac_2edit"]
selected_condition = "37_2hr"
time_enzyme = "lambdaN"
time_concentration = "250nM"
time_order = ["37_30min", "37_1hr", "37_2hr"]
# One processing pass serves both figures; each has its own correction family.
panel_definitions = [
    ("recorder_edit_summary", ["tada_type", "tada_conc", "edit_type"], 6),
    ("supp_recorder_time", ["tada_type", "tada_conc", "condition", "edit_type"], 3),
]
excluded_layout = "spacer5_0"
excluded_window = "1_3"
layout_pattern = r"^(spacer[35]_(?:0|10))"
expected_pairs = 15
edit_columns = [f"num_{number}_C" for number in range(1, 8)]
statistic_columns = [
    "mut_mean_percent", "mut_se_percent", "wt_mean_percent", "wt_se_percent",
    "difference_percentage_points", "fold_change",
]

samples = pd.read_csv(sample_file)
concentration_samples = (
    samples["tada_type"].isin(enzyme_order)
    & samples["tada_conc"].isin(concentration_order)
    & samples["condition"].eq(selected_condition)
)
time_samples = (
    samples["tada_type"].eq(time_enzyme)
    & samples["tada_conc"].eq(time_concentration)
    & samples["condition"].isin(time_order)
)
if concentration_samples.sum() != 6 or time_samples.sum() != len(time_order):
    raise ValueError("Expected six concentration samples and three time-course samples")
samples = samples.loc[
    concentration_samples | time_samples,
    ["sample_id", "tada_type", "tada_conc", "condition"],
]
if (samples["sample_id"].duplicated().any()
        or samples.duplicated(["tada_type", "tada_conc", "condition"]).any()):
    raise ValueError("Expected exactly one sample per enzyme, concentration, and time")

annotations = pd.read_csv(annotation_file)
annotations = annotations.loc[
    annotations["variable_type"].eq("boxb") & annotations["g_depleted"].eq("no"),
    ["reverse_complement", "variable_subpos", "oligo_name"],
].rename(columns={"reverse_complement": "barcode"})
annotations["layout"] = annotations["oligo_name"].str.extract(layout_pattern)
if annotations["layout"].isna().any():
    raise ValueError("Unrecognized recorder layout in barcode annotations")
stems = pd.read_csv(stem_file)

# Read each selected compressed sample once and discard unrelated constructs early.
frames = []
for sample_id in samples["sample_id"]:
    counts = pd.read_csv(
        f"{input_directory}/{sample_id}.csv.gz",
        usecols=["sample_id", "barcode", "insert", "umi_counts"] + edit_columns,
    )
    if not counts["sample_id"].eq(sample_id).all():
        raise ValueError(f"Sample identifiers do not match {sample_id}")
    counts = counts.merge(annotations, on="barcode", validate="many_to_one")
    counts = counts.merge(stems, on=["variable_subpos", "insert"], validate="many_to_one")
    frames.append(counts)
    print(f"Read {sample_id}: {len(counts)} WT/MUT constructs", flush=True)

constructs = pd.concat(frames, ignore_index=True).merge(
    samples, on="sample_id", validate="many_to_one"
)
excluded = constructs["layout"].eq(excluded_layout) & constructs["variable_subpos"].eq(excluded_window)
excluded_counts = constructs.loc[excluded].groupby(["sample_id", "insert_type"]).size()
if len(excluded_counts) != len(samples) * len(insert_order) or not excluded_counts.eq(1).all():
    raise ValueError("Expected one failed WT construct and its matched MUT per sample")
constructs = constructs.loc[~excluded].copy()
counts = constructs[["umi_counts"] + edit_columns]
if (not np.isfinite(counts.to_numpy()).all() or counts.lt(0).any().any()
        or constructs["umi_counts"].le(0).any()
        or constructs[edit_columns].sum(axis=1).gt(constructs["umi_counts"]).any()):
    raise ValueError("Invalid UMI counts or editing counts")
constructs["pair_id"] = constructs["layout"] + "__" + constructs["variable_subpos"]
constructs["frac_1edit"] = constructs[edit_columns[0]] / constructs["umi_counts"]
constructs["frac_2edit"] = constructs[edit_columns[1:]].sum(axis=1) / constructs["umi_counts"]
constructs = constructs.melt(
    id_vars=["sample_id", "barcode", "insert", "variable_subpos", "layout", "pair_id",
             "tada_type", "tada_conc", "condition", "insert_type", "umi_counts"] + edit_columns,
    value_vars=edit_order, var_name="edit_type", value_name="fraction_edited",
)
for column, order in zip(
    ["tada_type", "tada_conc", "edit_type", "insert_type"],
    [enzyme_order, concentration_order, edit_order, insert_order],
):
    constructs[column] = pd.Categorical(constructs[column], categories=order, ordered=True)
if constructs.isna().any().any():
    raise ValueError("Missing or unrecognized construct annotations")
construct_counts = constructs.groupby(["sample_id", "edit_type", "insert_type"], observed=True).size()
if (len(construct_counts) != len(samples) * len(edit_order) * len(insert_order)
        or not construct_counts.eq(expected_pairs).all()):
    raise ValueError("Expected 15 constructs per sample, edit category, and stem type")
# Summarize each panel independently from the same annotated count data.
os.makedirs(output_directory, exist_ok=True)
for panel_name, group_columns, expected_samples in panel_definitions:
    output_prefix = f"{output_directory}/{panel_name}"
    if panel_name == "recorder_edit_summary":
        panel_constructs = constructs.loc[
            constructs["condition"].eq(selected_condition)
        ].drop(columns="condition").copy()
    else:
        panel_constructs = constructs.loc[
            constructs["tada_type"].eq(time_enzyme)
            & constructs["tada_conc"].eq(time_concentration)
            & constructs["condition"].isin(time_order)
        ].copy()
        panel_constructs["condition"] = pd.Categorical(
            panel_constructs["condition"], categories=time_order, ordered=True
        )
    if panel_constructs["sample_id"].nunique() != expected_samples:
        raise ValueError(f"Unexpected sample count for {panel_name}")
    paired = panel_constructs.pivot(
        index=group_columns + ["pair_id"], columns="insert_type", values="fraction_edited"
    )
    if paired.isna().any().any():
        raise ValueError("Every retained WT construct must have a matched MUT")

    summary = panel_constructs.groupby(group_columns + ["insert_type"], observed=True)["fraction_edited"].agg(
        mean="mean", sd="std", n="size"
    ).reset_index()
    summary["se"] = 100 * summary.pop("sd") / np.sqrt(summary["n"])
    summary["mean"] *= 100
    summary = summary[group_columns + ["n", "insert_type", "mean", "se"]]

    tests = []
    for condition, pairs in paired.groupby(level=group_columns, observed=True):
        if len(pairs) != expected_pairs:
            raise ValueError(f"Incomplete pairs for {condition}")
        tests.append(dict(zip(group_columns, condition), p_value=ttest_rel(pairs["wt"], pairs["mut"]).pvalue))
    tests = pd.DataFrame(tests)
    tests["p_adjusted_bh"] = false_discovery_control(tests["p_value"].to_numpy(), method="bh")
    tests["significance"] = np.select(
        [tests["p_adjusted_bh"].lt(0.001), tests["p_adjusted_bh"].lt(0.01), tests["p_adjusted_bh"].lt(0.05)],
        ["***", "**", "*"], default="ns",
    )
    statistics = summary.pivot(index=group_columns, columns="insert_type", values=["mean", "se", "n"])
    statistics.columns = [f"{stem}_{measure}" for measure, stem in statistics.columns]
    statistics = statistics.reset_index().merge(tests, on=group_columns, validate="one_to_one")
    statistics["n_pairs"] = statistics[["mut_n", "wt_n"]].min(axis=1).astype(int)
    statistics["difference_percentage_points"] = statistics["wt_mean"] - statistics["mut_mean"]
    statistics["fold_change"] = statistics["wt_mean"] / statistics["mut_mean"]
    statistics = statistics.rename(columns={
        f"{stem}_{measure}": f"{stem}_{measure}_percent"
        for stem in insert_order for measure in ["mean", "se"]
    })[group_columns + ["n_pairs"] + statistic_columns + ["p_value", "p_adjusted_bh", "significance"]]
    if (len(summary) != expected_samples * len(edit_order) * len(insert_order)
            or len(statistics) != expected_samples * len(edit_order)
            or not np.isfinite(statistics[statistic_columns].to_numpy()).all()):
        raise ValueError(f"Unexpected or non-finite summaries for {panel_name}")

    panel_constructs.sort_values(group_columns + ["pair_id", "insert_type"]).to_csv(f"{output_prefix}_constructs.csv", index=False)
    summary.to_csv(f"{output_prefix}_plot_data.csv", index=False)
    statistics_suffix = "statistics" if panel_name == "recorder_edit_summary" else "statistics_full"
    statistics.to_csv(f"{output_prefix}_{statistics_suffix}.csv", index=False)
    print(f"{panel_name}: wrote {len(panel_constructs)} measurements, {len(summary)} means, and {len(statistics)} comparisons")
