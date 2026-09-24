# BoxB in vitro sequencing

## Shared recorder presentation

- The recorder panels share [presentation settings](./scripts/recorder_plot_style.R). Both figures use pink TadA–λN WT measurements, green TadA8.20 WT measurements, and gray MUT measurements. All panel titles are black. Axis lines and tick marks are thin (0.25 mm) and gray.
- Figure 2B pools distances and recorder adenosines into editing-count categories. Concentrations are labelled across the top of the columns, with TadA–λN above TadA8.20 and each enzyme named once vertically along the right edge; each concentration has an independent y-scale shared between its two enzyme panels. The compact figure is 4.5 × 2.5 inches; light and dark shades represent 1 edit and 2+ edits. WT precedes MUT within each edit category, matching Figure 2C; x-axis tick marks are omitted.
- Figure 2C separates distances and adenosines. Light-to-dark enzyme-specific shades encode distance. Its point-only gray legend, titled “boxB-recorder distance,” lists 0, 10, 20, and 30 nt; its x-axis rows are labelled “boxB,” “A context,” and “A position.” WT and MUT are grouped closely within each context, with wider gaps between contexts.

## Figure 2B recorder

- Run [preparation](./scripts/prepare_figure2b_recorder.py) once when counts or annotations change, then run [plotting](./scripts/plot_figure2b_recorder.R) independently. Each plotting run produces only the recorder concentration panel, as PNG at 96 DPI and PDF.
- The barcode splitting, read validation, UMI deduplication, counting, and combination rules in [the sequencing workflow](./scripts/run_analysis.smk) are unchanged. Existing combined counts are sufficient; raw FASTQ processing does not need to be repeated.
- Preparation reads eight sample files once from [combined recorder counts](./data/summary_stats_combined/): six for Figure 2B and two additional time points for the supplement. Samples are selected through [sample metadata](./annotations/sample_info.csv), then joined to [barcode annotations](./annotations/barcode_annotations.csv) and [existing WT/MUT stem definitions](./tables/boxb_wt_mut_stems.csv). The stem definitions were produced by [shared preparation](./scripts/01_prepare_data.R); the recorder workflow does not read its large prepared data tables.
- Selection remains both enzymes at 125, 250, and 500 nM, incubated at 37 °C for two hours. The failed spacer5_0/window 1_3 WT construct and its matched MUT are excluded, leaving 15 complete pairs per condition.
- Bars are equally weighted means of construct editing fractions, with standard errors across constructs. Comparisons use two-sided paired t-tests and Benjamini–Hochberg correction across all 12 tests. Fold changes are ratios of WT and MUT means.
- Baseline compatibility: the category labelled “2+ edits” deliberately retains the original sum of counts for 2–7 edits. Eight-edit molecules remain omitted. Correcting this definition is a separate scientific change.
- The main Figure 2B plotting script produces only the concentration panel. The time-course supplement has its own plotting script. Superseded fixed-recorder scripts and exploratory plots have been retired. The randomized-recorder analyses remain separate.

## Supplementary recorder time course

- Run [the time-course plotting script](./scripts/plot_recorder_time.R) after the same preparation step. It produces one supplementary panel with 30-minute, 1-hour, and 2-hour measurements of TadA–λN at 250 nM and 37 °C, using a common vertical scale and the original blue/gray recorder colors and bar layout.
- Use the same 15 matched WT/MUT pairs, failed-pair exclusion, mean and standard error, and 1-edit/2–7-edit definitions as Figure 2B. Apply Benjamini–Hochberg correction separately across the six supplementary paired comparisons. This replaces the earlier time-course calculation, which included the failed pair and used unpaired tests.
- Outputs: [180 construct measurements](./tables/supp_recorder_time_constructs.csv), [12 full-precision summaries](./tables/supp_recorder_time_plot_data.csv), and [six comparisons](./tables/supp_recorder_time_statistics_full.csv).
- Figure artifacts: [PNG](https://github.com/rasilab/bakker_2025/blob/master/analysis/boxb_in_vitro_sequencing/figures/supp_recorder_time.png?raw=1) and [PDF](https://github.com/rasilab/bakker_2025/blob/master/analysis/boxb_in_vitro_sequencing/figures/supp_recorder_time.pdf?raw=1).

![Supplementary recorder time course](https://github.com/rasilab/bakker_2025/blob/master/analysis/boxb_in_vitro_sequencing/figures/supp_recorder_time.png?raw=1)

## Figure 2C recorder positions

- [Generate Figure 2C](./scripts/plot_figure2c_recorder_positions.R) from the repository root. This is the distance-resolved WT/MUT comparison at 250 nM, 37 °C, and two hours. It replaces the separate pooled-by-position panel and the provisional distance-plot filenames.
- TadA–λN is the top facet with pink WT points and error bars; TadA8.20 is the bottom facet with green WT points and error bars. MUT measurements use gray shades in both facets. Within each enzyme, light-to-dark shades distinguish boxB recorder distances of 0, 10, 20, and 30 nt. Each facet has an independent y-axis range starting at zero, with 8% headroom above its highest error bar, and both share one bottom x-axis. There are no connecting lines.
- Each of the 128 points is an equally weighted construct mean ± SEM. The failed spacer5_0/window 1_3 pair is excluded for both enzymes, leaving three matched WT/MUT constructs at 0 nt and four at each other distance. These are construct measurements from one reaction per enzyme, not independent biological replicates.
- The script reads only the required columns of two compressed sample-count files, joins [sample metadata](./annotations/sample_info.csv), [barcode annotations](./annotations/barcode_annotations.csv), and [WT/MUT stem definitions](./tables/boxb_wt_mut_stems.csv), and uses [recorder-position annotations](./annotations/recorder_positions.csv). The position/context and count-column mappings were preserved from the previous position plot. Neither the retired pooled summary nor the large shared prepared table is a plotting dependency.
- The [summary table](./tables/fig2c_recorder_positions.csv) contains all 128 full-precision means, SEMs, and sample sizes. Pooled significance stars are omitted because they do not test distance-specific groups.
- Validation: all 128 means and SEMs agree with independent calculations from the shared prepared table within 1e-12. Renaming the distance plot to Figure 2C leaves its numerical table byte-identical. The initial complete distance-script run took 10.66 seconds in the preferred node-local container.
- Figure artifacts: [PNG](https://github.com/rasilab/bakker_2025/blob/master/analysis/boxb_in_vitro_sequencing/figures/fig2c_recorder_positions.png?raw=1) and [PDF](https://github.com/rasilab/bakker_2025/blob/master/analysis/boxb_in_vitro_sequencing/figures/fig2c_recorder_positions.pdf?raw=1).

![Figure 2C recorder positions by distance](https://github.com/rasilab/bakker_2025/blob/master/analysis/boxb_in_vitro_sequencing/figures/fig2c_recorder_positions.png?raw=1)

## Preserved randomized-recorder analyses

- The older [distance analysis](./scripts/plot_randomized_recorder_distance_legacy.R) and [randomized-context heatmap analysis](./scripts/plot_randomized_recorder_context_legacy.R) were extracted from the retired multi-panel script without changing their calculations. They retain their legacy output names and read the shared prepared table; they are reserved for a separate analysis pass.
- The [legacy R workflow](./scripts/generate_figures.R), [notebook](./scripts/generate_figures.ipynb), shared preparation, randomized-recorder distance tables, heatmap figures, and spacer-scan drafts are retained. Only obsolete Figure 2B concentration blocks were removed from the R workflow and notebook.
- Removed superseded fixed-recorder layout, paired-10-nt, and per-position layout drafts from the active folder after making local backups. The current entry points are the Figure 2B preparation/plotting scripts and the Figure 2C plotting script above.

## Run from the repository root

- Figure 2B and the time-course supplement:

- Inside the analysis container:

  ```bash
  python analysis/boxb_in_vitro_sequencing/scripts/prepare_figure2b_recorder.py
  Rscript analysis/boxb_in_vitro_sequencing/scripts/plot_figure2b_recorder.R
  Rscript analysis/boxb_in_vitro_sequencing/scripts/plot_recorder_time.R
  ```

- Outside a container, use the repository-preferred node-local image:

  ```bash
  module load Singularity
  singularity exec -B /fh /var/tmp/rasi-spacer-scan.NpIkvx/r_python_2.4.2.sif \
    python analysis/boxb_in_vitro_sequencing/scripts/prepare_figure2b_recorder.py
  singularity exec -B /fh /var/tmp/rasi-spacer-scan.NpIkvx/r_python_2.4.2.sif \
    Rscript analysis/boxb_in_vitro_sequencing/scripts/plot_figure2b_recorder.R
  singularity exec -B /fh /var/tmp/rasi-spacer-scan.NpIkvx/r_python_2.4.2.sif \
    Rscript analysis/boxb_in_vitro_sequencing/scripts/plot_recorder_time.R
  ```

- Figure 2C, inside the analysis container:

  ```bash
  Rscript analysis/boxb_in_vitro_sequencing/scripts/plot_figure2c_recorder_positions.R
  ```

## Outputs

- [Construct measurements](./tables/fig2b_recorder_constructs.csv): 360 rows, retaining sample, barcode, insert, pair identity, UMI counts, edit counts, and editing fraction.
- [Plotting summaries](./tables/fig2b_recorder_plot_data.csv): 24 means and standard errors at full precision.
- [Full statistics](./tables/fig2b_recorder_statistics_full.csv): 12 comparisons, including fold changes, raw and adjusted p-values, and significance labels. This and the plotting summary are the only plotting inputs, approximately 4 KB combined.
- [Legacy summary](./tables/fig2b_recorder.csv) and [legacy statistics](./tables/fig2b_recorder_statistics.csv): existing column names and two-decimal rounding preserved.
- Figure artifacts: [PNG](https://github.com/rasilab/bakker_2025/blob/master/analysis/boxb_in_vitro_sequencing/figures/fig2b_recorder.png?raw=1) and [PDF](https://github.com/rasilab/bakker_2025/blob/master/analysis/boxb_in_vitro_sequencing/figures/fig2b_recorder.pdf?raw=1).

![Figure 2B recorder](https://github.com/rasilab/bakker_2025/blob/master/analysis/boxb_in_vitro_sequencing/figures/fig2b_recorder.png?raw=1)

## Validation

- Final cleanup check: both regenerated PNGs are byte-identical to the approved figures, all Figure 2B tables are unchanged, and the Figure 2C table is byte-identical. The legacy randomized-context calculations and heatmap plotting blocks are unchanged, and all six active/retained R sources parse. Final container runs took 7.50 seconds for Figure 2B plotting and 8.66 seconds for Figure 2C preparation plus plotting.

- Compared all 360 construct measurements, 24 means and standard errors, and 12 comparisons against the original R implementation using the real data; agreement within absolute and relative tolerances of 1e-12. Existing rounded summaries also agree.
- Before the requested label update, the PNG and the PDF rendered at 144 DPI were pixel-identical to the saved baseline. The untethered enzyme is now labelled TadA8.20; numerical results are unchanged.
- Verified Figure 2B plotting from its small CSV inputs without sequencing summaries or the shared prepared table; the plotting script also sources the shared presentation settings.
- Independently recalculated the time-course means, standard errors, and paired tests from the three raw sample-count CSV files in R; all 12 summaries and six BH-corrected comparisons agree within 1e-12. Adding the supplement left all five Figure 2B CSV files byte-identical.
- Before replacing the pooled position panel, regenerated it with repository-root-relative paths: all 32 summary rows are unchanged and report n = 15; the PNG and PDF rendered at 144 DPI are pixel-identical to the recent WT/MUT plot. The legacy R sources parse after removal of the obsolete plotting blocks.
- Initial six-sample refactor elapsed times, including container startup: preparation 6.76 seconds; plotting 7.49 seconds. Peak memory was approximately 155 MiB and 200 MiB, respectively. These predate the two added time-course samples and are not runtime guarantees.

## Changelog

- **2026-09-24:** Finalized Figure 2B/2C presentation and the distance-resolved Figure 2C workflow; preserved numerical results, consolidated plotting styles and position annotations, retired superseded fixed-recorder code and figures, and retained randomized-recorder distance/heatmap analyses in separate legacy scripts.
- **2026-09-24:** Registered the WT/MUT recorder-position comparison as Figure 2C, documented 15 construct measurements per error bar, and retired the old single-enzyme plot and its plotting code while preserving other analyses.
- **2026-09-24:** Added a separately plotted TadA–λN time-course supplement; shared the preparation pass with Figure 2B while keeping its five CSV outputs byte-identical. The supplementary comparison uses the same filtering and paired-test method as Figure 2B.
- **2026-09-24:** Separated Figure 2B recorder preparation into Python and plotting into R; preserved numerical results and plot appearance apart from the requested TadA8.20 label; removed the bundled time course and its output, and cleaned up redundant filtered-recorder files locally. Other figure analyses and raw-read processing remain unchanged.
