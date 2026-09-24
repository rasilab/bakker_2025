"""Prepare the unchanged Figure 2B recorder comparison; run from the repository root."""

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
output_prefix = f"{output_directory}/fig2b_recorder"
enzyme_order = ["lambdaN", "tada_only"]
concentration_order = ["125nM", "250nM", "500nM"]
insert_order = ["mut", "wt"]
edit_order = ["frac_1edit", "frac_2edit"]
selected_condition = "37_2hr"
excluded_layout = "spacer5_0"
excluded_window = "1_3"
layout_pattern = r"^(spacer[35]_(?:0|10))"
expected_pairs = 15
edit_columns = [f"num_{number}_C" for number in range(1, 8)]
group_columns = ["tada_type", "tada_conc", "edit_type"]
round_columns = [
    "mut_mean_percent", "mut_se_percent", "wt_mean_percent", "wt_se_percent",
    "difference_percentage_points", "fold_change",
]

samples = pd.read_csv(sample_file)
samples = samples.loc[
    samples["tada_type"].isin(enzyme_order)
    & samples["tada_conc"].isin(concentration_order)
    & samples["condition"].eq(selected_condition),
    ["sample_id", "tada_type", "tada_conc"],
]
if (len(samples) != len(enzyme_order) * len(concentration_order)
        or samples["sample_id"].duplicated().any()
        or samples.duplicated(["tada_type", "tada_conc"]).any()):
    raise ValueError("Expected exactly one sample for each of the six conditions")

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
             "tada_type", "tada_conc", "insert_type", "umi_counts"] + edit_columns,
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
paired = constructs.pivot(
    index=group_columns + ["pair_id"], columns="insert_type", values="fraction_edited"
)
if paired.isna().any().any():
    raise ValueError("Every retained WT construct must have a matched MUT")

summary = constructs.groupby(group_columns + ["insert_type"], observed=True)["fraction_edited"].agg(
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
})[group_columns + ["n_pairs"] + round_columns + ["p_value", "p_adjusted_bh", "significance"]]
if len(summary) != 24 or len(statistics) != 12 or not np.isfinite(statistics[round_columns].to_numpy()).all():
    raise ValueError("Expected 24 finite means and 12 finite comparisons")

os.makedirs(output_directory, exist_ok=True)
constructs.sort_values(group_columns + ["pair_id", "insert_type"]).to_csv(f"{output_prefix}_constructs.csv", index=False)
summary.to_csv(f"{output_prefix}_plot_data.csv", index=False)
statistics.to_csv(f"{output_prefix}_statistics_full.csv", index=False)
# Retain the existing public summary schemas and rounding for downstream users.
summary.round({"mean": 2, "se": 2}).to_csv(f"{output_prefix}.csv", index=False)
statistics.round(dict.fromkeys(round_columns, 2)).to_csv(f"{output_prefix}_statistics.csv", index=False)
print(f"Wrote {len(constructs)} measurements, {len(summary)} means, and {len(statistics)} comparisons")
