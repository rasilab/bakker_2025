"""Plot Figure 2B recorder editing separately by reporter layout."""

import os

os.environ.setdefault("MPLCONFIGDIR", os.path.expanduser("~/.cache/matplotlib"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

# Analysis choices
input_directory = "analysis/boxb_in_vitro_sequencing/data/summary_stats_combined"
annotation_file = "analysis/boxb_in_vitro_sequencing/annotations/barcode_annotations.csv"
sample_info_file = "analysis/boxb_in_vitro_sequencing/annotations/sample_info.csv"
stem_file = "analysis/boxb_in_vitro_sequencing/tables/boxb_wt_mut_stems.csv"
output_png = "analysis/boxb_in_vitro_sequencing/figures/fig2b_by_layout.png"
output_pdf = "analysis/boxb_in_vitro_sequencing/figures/fig2b_by_layout.pdf"
figure_width = 4.7
figure_height = 4.0
font_family = "Helvetica"
font_size = 7
enzyme_order = ["tada_only", "lambdaN"]
enzyme_labels = {"tada_only": "TadA", "lambdaN": "λN–TadA"}
concentration_order = ["125nM", "250nM", "500nM"]
concentration_labels = {"125nM": "125 nM", "250nM": "250 nM", "500nM": "500 nM"}
edit_order = ["1_edit", "2plus_edits"]
edit_labels = {"1_edit": "1 edit", "2plus_edits": "2+ edits"}
insert_order = ["wt", "mut"]
insert_labels = {"wt": "WT", "mut": "MUT"}
insert_offsets = {"wt": -0.16, "mut": 0.16}
insert_markers = {"wt": "o", "mut": "^"}
layout_order = ["spacer5_0", "spacer5_10", "spacer3_10", "spacer3_0"]
layout_labels = {
    "spacer5_0": "0 nt",
    "spacer5_10": "10 nt",
    "spacer3_10": "20 nt",
    "spacer3_0": "30 nt",
}
layout_offsets = dict(zip(layout_order, [-0.045, -0.015, 0.015, 0.045]))
layout_colors = dict(zip(layout_order, ["#0072B2", "#E69F00", "#009E73", "#CC79A7"]))

plt.rcParams.update({
    "font.family": font_family,
    "font.size": font_size,
    "axes.labelsize": font_size,
    "xtick.labelsize": font_size,
    "ytick.labelsize": font_size,
    "legend.fontsize": font_size,
    "axes.linewidth": 0.5,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
    "xtick.major.size": 2,
    "ytick.major.size": 2,
    "lines.linewidth": 0.5,
    "pdf.fonttype": 42,
})

annotations = pd.read_csv(annotation_file)
annotations = annotations.loc[
    annotations["variable_type"].eq("boxb") & annotations["g_depleted"].eq("no"),
    ["reverse_complement", "variable_subpos", "oligo_name"],
].rename(columns={"reverse_complement": "barcode"})
annotations["layout"] = annotations["oligo_name"].str.extract(r"^(spacer[35]_(?:0|10))")
stems = pd.read_csv(stem_file)
samples = pd.read_csv(sample_info_file)
samples = samples.loc[
    samples["tada_type"].isin(enzyme_order)
    & samples["tada_conc"].isin(concentration_order)
    & samples["condition"].eq("37_2hr"),
    ["sample_id", "tada_type", "tada_conc"],
]

target_columns = ["sample_id", "barcode", "insert", "umi_counts"] + [
    f"num_{edit}_C" for edit in range(1, 9)
]
sample_frames = []
for sample_id in samples["sample_id"]:
    sample_frames.append(
        pd.read_csv(f"{input_directory}/{sample_id}.csv.gz", usecols=target_columns)
    )
plot_data = pd.concat(sample_frames, ignore_index=True)
plot_data = plot_data.merge(annotations, on="barcode", how="inner")
plot_data = plot_data.merge(stems, on=["variable_subpos", "insert"], how="inner")
plot_data = plot_data.merge(samples, on="sample_id", how="inner")
plot_data["fraction_1_edit"] = plot_data["num_1_C"] / plot_data["umi_counts"]
plot_data["fraction_2plus_edits"] = (
    plot_data[[f"num_{edit}_C" for edit in range(2, 9)]].sum(axis=1)
    / plot_data["umi_counts"]
)
plot_data = plot_data.melt(
    id_vars=["sample_id", "tada_type", "tada_conc", "layout", "insert_type"],
    value_vars=["fraction_1_edit", "fraction_2plus_edits"],
    var_name="edit_type",
    value_name="fraction_edited",
)
plot_data["edit_type"] = plot_data["edit_type"].str.removeprefix("fraction_")

summary = (
    plot_data.groupby(
        ["sample_id", "tada_type", "tada_conc", "edit_type", "insert_type", "layout"],
        observed=True,
    )["fraction_edited"]
    .agg(mean="mean", se=lambda values: values.std(ddof=1) / np.sqrt(len(values)), n="size")
    .reset_index()
)
if not summary["n"].eq(4).all():
    raise ValueError("Expected four stem-window measurements per region, layout, and insert type")

fig, axes = plt.subplots(
    len(concentration_order),
    len(enzyme_order),
    figsize=(figure_width, figure_height),
    sharex=True,
    sharey="row",
    layout="constrained",
)
x_positions = np.arange(len(edit_order))
insert_tick_positions = [
    edit_index + insert_offsets[insert_type]
    for edit_index in range(len(edit_order))
    for insert_type in insert_order
]
insert_tick_labels = [
    insert_labels[insert_type]
    for _ in edit_order
    for insert_type in insert_order
]

for row, concentration in enumerate(concentration_order):
    for column, enzyme in enumerate(enzyme_order):
        axis = axes[row, column]
        facet_data = summary.loc[
            summary["tada_conc"].eq(concentration) & summary["tada_type"].eq(enzyme)
        ]
        for edit_index, edit_type in enumerate(edit_order):
            for insert_type in insert_order:
                group_data = (
                    facet_data.loc[
                        facet_data["edit_type"].eq(edit_type)
                        & facet_data["insert_type"].eq(insert_type)
                    ]
                    .set_index("layout")
                    .loc[layout_order]
                )
                plot_x = np.array([
                    edit_index + insert_offsets[insert_type] + layout_offsets[layout]
                    for layout in layout_order
                ])
                for layout, point_x, mean, se in zip(
                    layout_order, plot_x, group_data["mean"], group_data["se"]
                ):
                    axis.errorbar(
                        point_x,
                        mean * 100,
                        yerr=se * 100,
                        color=layout_colors[layout],
                        marker=insert_markers[insert_type],
                        markersize=2.5,
                        markeredgewidth=0,
                        linewidth=0,
                        elinewidth=0.4,
                        capsize=1,
                        zorder=3,
                    )
        axis.set_xlim(-0.42, 1.42)
        axis.set_xticks(insert_tick_positions, insert_tick_labels)
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
        if row == 0:
            axis.text(0.5, 1.03, enzyme_labels[enzyme], transform=axis.transAxes,
                      ha="center", fontweight="bold")
        if column == 0:
            axis.text(-0.28, 0.5, concentration_labels[concentration], transform=axis.transAxes,
                      ha="center", va="center", rotation=90, fontweight="bold")
        if row < len(concentration_order) - 1:
            axis.tick_params(labelbottom=False)
        else:
            for edit_index, edit_type in enumerate(edit_order):
                axis.text(
                    edit_index,
                    -0.19,
                    edit_labels[edit_type],
                    transform=axis.get_xaxis_transform(),
                    ha="center",
                    va="top",
                )

fig.supxlabel("Recorder editing category")
fig.supylabel("Edited RNA (%)")
legend_handles = [
    Line2D([], [], color=layout_colors[layout], marker="o", linestyle="none",
           markersize=3, markeredgewidth=0, label=layout_labels[layout])
    for layout in layout_order
]
fig.legend(handles=legend_handles, loc="outside right center", frameon=False, title="Distance")
fig.savefig(output_png, format="png", dpi=96, bbox_inches="tight", pad_inches=0.05, transparent=True)
fig.savefig(output_pdf, format="pdf", bbox_inches="tight", pad_inches=0.05, transparent=True)
print(summary.to_string(index=False))
print(f"Wrote {output_png} and {output_pdf}")
