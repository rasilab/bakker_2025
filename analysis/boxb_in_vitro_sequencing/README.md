# BoxB in vitro sequencing

- Three main recorder plots use descriptive filenames; panel letters are assigned during manuscript assembly.
- [Shared colors and typography](./scripts/recorder_plot_style.R): pink TadA–λN, green TadA8.20, gray MUT, Helvetica. Standalone figures retain PNG at 96 DPI and PDF; the heat map also retains SVG.
- Raw counts, sequencing workflows, Figure 3/4 loop and stem plots, recruitment analyses, and later figures are preserved. The obsolete aggregate Figure 2 loop plot was already removed on September 23.

## Recorder edit summary

- [Prepare](./scripts/prepare_recorder_edit_summary.py), then [plot](./scripts/plot_recorder_edit_summary.R). Preparation reads eight combined-count samples once and also supplies the time-course supplement.
- Both enzymes at 125, 250, and 500 nM, 37 °C, two hours. Exclude the failed spacer5_0/window 1_3 WT/MUT pair; retain 15 matched construct pairs per condition.
- Bars show equally weighted construct means ± SEM. Tests are paired, two-sided t-tests, with Benjamini–Hochberg correction across 12 comparisons. The existing “2+ edits” category remains the sum of 2–7 edits; eight-edit molecules remain omitted.
- Tables: [360 construct measurements](./tables/recorder_edit_summary_constructs.csv), [24 means and SEMs](./tables/recorder_edit_summary_plot_data.csv), [12 full-precision comparisons](./tables/recorder_edit_summary_statistics.csv). Redundant rounded exports are retired.

## Recorder sites and distances

- [Entry script](./scripts/plot_recorder_site_distance.R) uses [shared preparation and plotting](./scripts/recorder_site_distance_panel.R) for the main 250 nM plot and separate 125/500 nM supplements.
- Eight annotated adenosines, four recorder–boxB distances, WT/MUT and both enzymes at 37 °C for two hours. [The table](./tables/recorder_site_distance.csv) retains 128 means, SEMs and sample sizes.
- The failed pair is excluded: three matched constructs at 0 nt, four at each other distance. Measurements are construct variation within one reaction per enzyme, not independent reaction replicates. Enzyme-specific color shades encode 0, 10, 20 and 30 nt.

## Recorder context heat map

- [Prepare](./scripts/prepare_recorder_context_heatmap.py), then [plot](./scripts/plot_recorder_context_heatmap.R), using [sample/site mappings](./annotations/recorder_context_selection.csv) and [selected geometries](./annotations/recorder_context_heatmap_selection.csv).
- TadA–λN and TadA8.20 at 250 nM, 37 °C, two hours. Pool A2 (UAG), A8 (UAC), A10 (CAC), and A13 (CAU), with recorder 5′ of boxB at 0/10 nt and 3′ at 20/30 nt. Preparation verifies the same geometry set and 17-nt recorder length as the boxB-variation library; barcodes differ between libraries.
- Each tile is 100 × summed edited counts / summed site-specific UMI coverage across four sites, four geometries and 64 sequence backgrounds (1,024 site/variant observations). This is a percentage of edited site observations; shared molecules contribute once per evaluated site. Observed flanking identities remain separate, without an editing correction or coverage cutoff.
- Tables retain [512 site/configuration measurements](./tables/recorder_context_heatmap_sites.csv) and [32 pooled contexts](./tables/recorder_context_heatmap.csv), with full-precision fractions and count totals.
- Two side-by-side heat maps retain the shared 0–40% scale, integer percentage labels and no color bars. Rows give the 5′ flank and columns the 3′ flank; dimensions are 3.2 × 1.8 inches.
- UAU remains highest: 32.50% versus 20.78% for UAA with TadA–λN, and 34.56% versus 21.80% with TadA8.20. This ranking persists with equal site/configuration weighting and after omitting any one site or geometry. These are not independent biological reaction replicates.
- Archived validation also checked 3′ spacers at 10/20/30 nt: UAU led all six pooled enzyme/geometry comparisons and 22 of 24 individual site comparisons. At A8 with a 10-nt spacer, UAC narrowly led for both enzymes.

## Retained supplements

- [Time-course plotting](./scripts/plot_recorder_time.R) uses the edit-summary preparation: TadA–λN at 250 nM, 37 °C, 30 minutes/one hour/two hours; the same 15 pairs and edit categories. Its six paired tests have a separate Benjamini–Hochberg correction family.
- Time-course tables: [construct measurements](./tables/supp_recorder_time_constructs.csv), [plot data](./tables/supp_recorder_time_plot_data.csv), [statistics](./tables/supp_recorder_time_statistics_full.csv).
- Site-distance supplements use the same entry script and method as the main plot: [125 nM table](./tables/supp_recorder_positions_125nM.csv), [500 nM table](./tables/supp_recorder_positions_500nM.csv).

## Run from the repository root

- Run inside the analysis container, using the node-local image specified in [repository instructions](../../AGENTS.md):

  ```bash
  python analysis/boxb_in_vitro_sequencing/scripts/prepare_recorder_edit_summary.py
  Rscript analysis/boxb_in_vitro_sequencing/scripts/plot_recorder_edit_summary.R
  Rscript analysis/boxb_in_vitro_sequencing/scripts/plot_recorder_time.R
  Rscript analysis/boxb_in_vitro_sequencing/scripts/plot_recorder_site_distance.R
  python analysis/boxb_in_vitro_sequencing/scripts/prepare_recorder_context_heatmap.py
  Rscript analysis/boxb_in_vitro_sequencing/scripts/plot_recorder_context_heatmap.R
  ```

## Figure artifacts

- Recorder edit summary: [PNG](https://github.com/rasilab/bakker_2025/blob/master/analysis/boxb_in_vitro_sequencing/figures/recorder_edit_summary.png?raw=1), [PDF](https://github.com/rasilab/bakker_2025/blob/master/analysis/boxb_in_vitro_sequencing/figures/recorder_edit_summary.pdf?raw=1).

![Recorder edit summary](https://github.com/rasilab/bakker_2025/blob/master/analysis/boxb_in_vitro_sequencing/figures/recorder_edit_summary.png?raw=1)

- Recorder sites and distances: [PNG](https://github.com/rasilab/bakker_2025/blob/master/analysis/boxb_in_vitro_sequencing/figures/recorder_site_distance.png?raw=1), [PDF](https://github.com/rasilab/bakker_2025/blob/master/analysis/boxb_in_vitro_sequencing/figures/recorder_site_distance.pdf?raw=1).

![Recorder sites and distances](https://github.com/rasilab/bakker_2025/blob/master/analysis/boxb_in_vitro_sequencing/figures/recorder_site_distance.png?raw=1)

- Recorder context heat map: [PNG](https://github.com/rasilab/bakker_2025/blob/master/analysis/boxb_in_vitro_sequencing/figures/recorder_context_heatmap.png?raw=1), [PDF](https://github.com/rasilab/bakker_2025/blob/master/analysis/boxb_in_vitro_sequencing/figures/recorder_context_heatmap.pdf?raw=1), [SVG](https://github.com/rasilab/bakker_2025/blob/master/analysis/boxb_in_vitro_sequencing/figures/recorder_context_heatmap.svg?raw=1).

![Recorder context heat map](https://github.com/rasilab/bakker_2025/blob/master/analysis/boxb_in_vitro_sequencing/figures/recorder_context_heatmap.png?raw=1)

- Recorder time course: [PNG](https://github.com/rasilab/bakker_2025/blob/master/analysis/boxb_in_vitro_sequencing/figures/supp_recorder_time.png?raw=1), [PDF](https://github.com/rasilab/bakker_2025/blob/master/analysis/boxb_in_vitro_sequencing/figures/supp_recorder_time.pdf?raw=1).

![Recorder time course](https://github.com/rasilab/bakker_2025/blob/master/analysis/boxb_in_vitro_sequencing/figures/supp_recorder_time.png?raw=1)

- Recorder sites at 125 nM: [PNG](https://github.com/rasilab/bakker_2025/blob/master/analysis/boxb_in_vitro_sequencing/figures/supp_recorder_positions_125nM.png?raw=1), [PDF](https://github.com/rasilab/bakker_2025/blob/master/analysis/boxb_in_vitro_sequencing/figures/supp_recorder_positions_125nM.pdf?raw=1).

![Recorder sites at 125 nM](https://github.com/rasilab/bakker_2025/blob/master/analysis/boxb_in_vitro_sequencing/figures/supp_recorder_positions_125nM.png?raw=1)

- Recorder sites at 500 nM: [PNG](https://github.com/rasilab/bakker_2025/blob/master/analysis/boxb_in_vitro_sequencing/figures/supp_recorder_positions_500nM.png?raw=1), [PDF](https://github.com/rasilab/bakker_2025/blob/master/analysis/boxb_in_vitro_sequencing/figures/supp_recorder_positions_500nM.pdf?raw=1).

![Recorder sites at 500 nM](https://github.com/rasilab/bakker_2025/blob/master/analysis/boxb_in_vitro_sequencing/figures/supp_recorder_positions_500nM.png?raw=1)

## Validation and cleanup

- Cleanup verification: all 11 numerical tables and the selection CSV are byte-identical; all six PNGs are byte-identical and all six PDFs render pixel-identically at 144 DPI. All 14 retained R sources and every notebook code cell parse.
- Previous independent source-count calculations reproduced the recorder summaries and statistics within 1e-12 and all 32 heat-map numerators and denominators exactly.
- [Cleanup audit](../../specs/recorder_figure_cleanup_audit.md) and [file inventory](../../specs/recorder_figure_cleanup_inventory.csv) document renamed and retired artifacts. Uncommitted explorations and geometry-validation dependencies were archived before removal.
- The mixed R workflow and notebook retain subsequent-figure calculations, including shared constant-context and single-flank supplementary calculations.

## Changelog

- **2026-09-26:** Applied the approved recorder cleanup: descriptive names for three main plots, removed obsolete scripts/tables/figures and redundant exports, archived exploratory validation, and preserved recorder supplements and Figure 3/4 loop plots.

- **2026-09-25:** Removed redundant color bars from the pooled-context heat map and reduced its width to 3.2 inches. Integer labels show percent edited; the shared 0–40% scale and full-precision data are unchanged.

- **2026-09-25:** Updated Figure 2E to pool the same four recorder geometries as the September 24 Figure 2B/2C workflow: 5′/0, 5′/10, 3′/20, and 3′/30 nucleotides. Added direct source-count preparation with recorder-length and geometry checks, retained 512 site/configuration measurements, independently verified all 32 pooled tiles, and retained integer labels and the shared pink/green palette. Expanded the separate geometry check to include 3′/30 nucleotides.

- **2026-09-25:** Checked Figure 2E context preference at recorder 3′ of boxB with 10- and 20-nucleotide spacers. UAU remains the top pooled context for both enzymes under count pooling, equal-site weighting, and omission of each site. Saved source-derived site/context and comparison tables, verified overlap with the earlier distance analysis, and documented A8 exceptions at the 10-nucleotide spacer.

- **2026-09-25:** Rounded Figure 2E tile labels to integer percentages, preserving full-precision values and colors. Checked UAU rankings at all eight enzyme/site combinations and under equal-site weighting and omission of each site; documented the near tie at A8 for TadA8.20 and the lack of independent reaction replication.

- **2026-09-25:** Adapted Figure 2E to compare 250 nM TadA–λN and TadA8.20 using the recent pooled-context measurements at one fixed geometry: 10-nucleotide spacer, recorder 5′ of boxB. Each of 16 contexts pools counts across A2/A8/A10/A13. Added the selection and full-precision plotting table, matched the Figure 2B/2C palette, and generated side-by-side heat maps in PNG/PDF/SVG.

- **2026-09-24:** Added the pooled context-distance scan for UAA/UAC/UAG/UAU across 0–30-nucleotide spacers, retaining separate 5′/3′ orientations and four site markers. The two enzyme rows share a y-axis scale and use pink/green orientation shades. All 1,024 pooled measurements have complete 64-variant coverage; independently checked source counts and agreement with the finalized fixed-distance plot.

- **2026-09-24:** Added direct preparation and plotting for the 250 nM TadA–λN/TadA8.20 context comparison at a fixed 10-nucleotide spacer, recorder 5′ of boxB. Each of 16 contexts shows pooled-UMI percentages for A2, A8, A10, and A13 plus their equally weighted mean bar. Added count audit tables, a shared “Edited RNA percentage” label, and enzyme panel titles; independently validated all 128 pooled points and all 32 mean bars.

- **2026-09-24:** Documented the legacy Figure 2E heat-map workflow with a reproducible SVG/PNG/PDF schematic, source-linked explanation, and an independent count-table audit. The analysis and existing heat-map outputs are unchanged.

- **2026-09-24:** Made the upper-facet zero reference lines thinner, light gray, and dotted in Figure 2C and both supplementary concentration plots. Pixel comparisons confirm that the lower panels are unchanged; all numerical tables remain byte-identical.

- **2026-09-24:** Extended the Figure 2C entry script to generate matching 125 nM and 500 nM supplementary panels through shared preparation/plotting code, using descriptive filenames until supplement numbering is assigned. Added a thin gray zero baseline to the upper facet of all three recorder-position figures.

- **2026-09-24:** Finalized Figure 2B/2C presentation and the distance-resolved Figure 2C workflow; preserved numerical results, consolidated plotting styles and position annotations, retired superseded fixed-recorder code and figures, and retained randomized-recorder distance/heatmap analyses in separate legacy scripts.
- **2026-09-24:** Registered the WT/MUT recorder-position comparison as Figure 2C, documented 15 construct measurements per error bar, and retired the old single-enzyme plot and its plotting code while preserving other analyses.
- **2026-09-24:** Added a separately plotted TadA–λN time-course supplement; shared the preparation pass with Figure 2B while keeping its five CSV outputs byte-identical. The supplementary comparison uses the same filtering and paired-test method as Figure 2B.
- **2026-09-24:** Separated Figure 2B recorder preparation into Python and plotting into R; preserved numerical results and plot appearance apart from the requested TadA8.20 label; removed the bundled time course and its output, and cleaned up redundant filtered-recorder files locally. Other figure analyses and raw-read processing remain unchanged.
