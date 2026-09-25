"""Plot four recorder-site percentages per context after pooling background UMIs."""

import os
os.environ.setdefault("MPLCONFIGDIR", os.path.expanduser("~/.cache/matplotlib"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

# Run from the repository root; plotting reads only the prepared summary CSV.
ROOT = "analysis/boxb_in_vitro_sequencing"
SUMMARY_FILE = f"{ROOT}/tables/recorder_context_summary.csv"
OUTPUT = f"{ROOT}/figures/recorder_context_variants"
FIGURE_SIZE = (6.2, 3.6)
FONT_FAMILY, FONT_SIZE = "Helvetica", 7
LEGEND_COLOR, AXIS_COLOR = "#777777", "#999999"
POINT_SIZE, POSITION_OFFSET, Y_HEADROOM = 12, 0.24, 1.08
BAR_WIDTH, BAR_ALPHA = 0.72, 0.3
RNA_BASES = "ACGU"
EXPECTED_SAMPLES, EXPECTED_POSITIONS, EXPECTED_VARIANTS = 2, 4, 64

plt.rcParams.update({
    "font.family": FONT_FAMILY, "font.size": FONT_SIZE,
    "axes.labelsize": FONT_SIZE, "axes.titlesize": FONT_SIZE,
    "xtick.labelsize": FONT_SIZE, "ytick.labelsize": FONT_SIZE,
    "legend.fontsize": FONT_SIZE, "legend.title_fontsize": FONT_SIZE,
    "figure.labelsize": FONT_SIZE,
    "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2, "ytick.major.size": 2, "lines.linewidth": 0.5,
    "pdf.fonttype": 42,
})
summary = pd.read_csv(SUMMARY_FILE)
identity_columns = ["sample_id", "target_pos_to_boxb", "target_dist", "tada_type", "tada_conc",
                    "condition", "enzyme_label", "panel_order", "point_color"]
identity = summary[identity_columns].drop_duplicates().sort_values("panel_order")
shared_columns = ["target_pos_to_boxb", "target_dist", "tada_conc", "condition"]
site_columns = ["position_id", "reference_context", "position_order", "marker"]
sites = summary[site_columns].drop_duplicates().sort_values("position_order")
if len(identity) != EXPECTED_SAMPLES or len(identity[shared_columns].drop_duplicates()) != 1:
    raise ValueError("Expected two enzyme samples with the same concentration, incubation, and geometry")
if len(sites) != EXPECTED_POSITIONS:
    raise ValueError("Expected four recorder-site definitions")
central_bases = sites.reference_context.str[1].unique()
if len(central_bases) != 1:
    raise ValueError("Selected sites do not share the same central base")
context_order = [f"{five}{central_bases[0]}{three}" for five in RNA_BASES for three in RNA_BASES]
keys = ["sample_id", "position_id", "rna_context"]
if (len(summary) != EXPECTED_SAMPLES * EXPECTED_POSITIONS * len(context_order)
        or summary.duplicated(keys).any() or not summary.n_variants.eq(EXPECTED_VARIANTS).all()):
    raise ValueError("Each enzyme/site/context must contain 64 measured variants")
if (not summary.total_umis.gt(0).all()
        or not summary.total_edited_umis.between(0, summary.total_umis).all()
        or not np.allclose(summary.pooled_fraction_central_a_edited,
                           summary.total_edited_umis / summary.total_umis, rtol=0, atol=1e-12)):
    raise ValueError("Plot fractions do not match pooled edited and total UMI counts")

fig, axes = plt.subplots(EXPECTED_SAMPLES, 1, figsize=FIGURE_SIZE,
                         layout="constrained", sharex=True, sharey=True)
x = np.arange(len(context_order))
offsets = np.linspace(-POSITION_OFFSET, POSITION_OFFSET, EXPECTED_POSITIONS)
for axis, choice in zip(axes, identity.itertuples()):
    site_values = summary.loc[summary.sample_id.eq(choice.sample_id)].pivot(
        index="rna_context", columns="position_id", values="pooled_fraction_central_a_edited")
    site_values = site_values.loc[context_order, sites.position_id]
    if site_values.isna().any().any():
        raise ValueError("Each mean bar requires all four pooled site percentages")
    axis.bar(x, site_values.mean(axis=1) * 100, width=BAR_WIDTH,
             color=choice.point_color, alpha=BAR_ALPHA, edgecolor="none", zorder=1)
    for offset, site in zip(offsets, sites.itertuples()):
        data = summary.loc[summary.sample_id.eq(choice.sample_id)
                           & summary.position_id.eq(site.position_id)].set_index("rna_context").loc[context_order]
        axis.scatter(x + offset, data.pooled_fraction_central_a_edited * 100,
                     s=POINT_SIZE, marker=site.marker, color=choice.point_color, edgecolors="none", zorder=3)
    axis.set_xticks(x, context_order)
    axis.set_xlim(-0.65, len(context_order) - 0.35)
    axis.set_ylim(0, summary.pooled_fraction_central_a_edited.max() * 100 * Y_HEADROOM)
    axis.set_title(choice.enzyme_label, color="black")
    axis.spines[["top", "right"]].set_visible(False)
    for spine in ["bottom", "left"]:
        axis.spines[spine].set_color(AXIS_COLOR)
    axis.tick_params(axis="x", length=0)
    axis.tick_params(axis="y", color=AXIS_COLOR)
axes[-1].set_xlabel("Observed RNA context (central A fixed)")
fig.supylabel("Edited RNA percentage")
fig.legend(handles=[Patch(facecolor=LEGEND_COLOR, alpha=BAR_ALPHA, label="Mean (4 sites)")]
                   + [Line2D([], [], color=LEGEND_COLOR, marker=site.marker, linestyle="none",
                           markersize=np.sqrt(POINT_SIZE), markeredgewidth=0, label=site.position_id)
                    for site in sites.itertuples()],
           loc="outside right center", frameon=False)
fig.savefig(f"{OUTPUT}.png", dpi=96, facecolor="white")
fig.savefig(f"{OUTPUT}.pdf", facecolor="white")
print(f"Wrote {OUTPUT}.png and {OUTPUT}.pdf ({len(summary)} site points; equal-weight four-site mean bars)")
print({extension: os.path.getsize(f"{OUTPUT}.{extension}") for extension in ["png", "pdf"]})
