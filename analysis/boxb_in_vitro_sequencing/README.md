# BoxB in vitro sequencing

## Figure 2B recorder

- Run [preparation](./scripts/prepare_figure2b_recorder.py) once when counts or annotations change, then run [plotting](./scripts/plot_figure2b_recorder.R) independently. Each plotting run produces only the recorder concentration panel, as PNG at 96 DPI and PDF.
- The barcode splitting, read validation, UMI deduplication, counting, and combination rules in [the sequencing workflow](./scripts/run_analysis.smk) are unchanged. Existing combined counts are sufficient; raw FASTQ processing does not need to be repeated.
- Preparation reads six sample files from [combined recorder counts](./data/summary_stats_combined/), selected through [sample metadata](./annotations/sample_info.csv), and joins [barcode annotations](./annotations/barcode_annotations.csv) and [existing WT/MUT stem definitions](./tables/boxb_wt_mut_stems.csv). The stem definitions were produced by [shared preparation](./scripts/01_prepare_data.R); the new Figure 2B workflow does not read its large prepared data tables.
- Selection remains both enzymes at 125, 250, and 500 nM, incubated at 37 °C for two hours. The failed spacer5_0/window 1_3 WT construct and its matched MUT are excluded, leaving 15 complete pairs per condition.
- Bars are equally weighted means of construct editing fractions, with standard errors across constructs. Comparisons use two-sided paired t-tests and Benjamini–Hochberg correction across all 12 tests. Fold changes are ratios of WT and MUT means.
- Baseline compatibility: the category labelled “2+ edits” deliberately retains the original sum of counts for 2–7 edits. Eight-edit molecules remain omitted. Correcting this definition is a separate scientific change.
- This Figure 2B workflow replaces the recorder-only script and produces no time-course figure. The redundant filtered-recorder script and its outputs were removed locally. Alternative recorder views and the older multi-panel workflow are outside this commit.

## Run from the repository root

- Inside the analysis container:

  ```bash
  python analysis/boxb_in_vitro_sequencing/scripts/prepare_figure2b_recorder.py
  Rscript analysis/boxb_in_vitro_sequencing/scripts/plot_figure2b_recorder.R
  ```

- Outside a container, use the repository-preferred node-local image:

  ```bash
  module load Singularity
  singularity exec -B /fh /var/tmp/rasi-spacer-scan.NpIkvx/r_python_2.4.2.sif \
    python analysis/boxb_in_vitro_sequencing/scripts/prepare_figure2b_recorder.py
  singularity exec -B /fh /var/tmp/rasi-spacer-scan.NpIkvx/r_python_2.4.2.sif \
    Rscript analysis/boxb_in_vitro_sequencing/scripts/plot_figure2b_recorder.R
  ```

## Outputs

- [Construct measurements](./tables/fig2b_recorder_constructs.csv): 360 rows, retaining sample, barcode, insert, pair identity, UMI counts, edit counts, and editing fraction.
- [Plotting summaries](./tables/fig2b_recorder_plot_data.csv): 24 means and standard errors at full precision.
- [Full statistics](./tables/fig2b_recorder_statistics_full.csv): 12 comparisons, including fold changes, raw and adjusted p-values, and significance labels. This and the plotting summary are the only plotting inputs, approximately 4 KB combined.
- [Legacy summary](./tables/fig2b_recorder.csv) and [legacy statistics](./tables/fig2b_recorder_statistics.csv): existing column names and two-decimal rounding preserved.
- Figure artifacts: [PNG](https://github.com/rasilab/bakker_2025/blob/master/analysis/boxb_in_vitro_sequencing/figures/fig2b_recorder.png?raw=1) and [PDF](https://github.com/rasilab/bakker_2025/blob/master/analysis/boxb_in_vitro_sequencing/figures/fig2b_recorder.pdf?raw=1).

![Figure 2B recorder](https://github.com/rasilab/bakker_2025/blob/master/analysis/boxb_in_vitro_sequencing/figures/fig2b_recorder.png?raw=1)

## Validation

- Compared all 360 construct measurements, 24 means and standard errors, and 12 comparisons against the original R implementation using the real data; agreement within absolute and relative tolerances of 1e-12. Existing rounded summaries also agree.
- Before the requested label update, the PNG and the PDF rendered at 144 DPI were pixel-identical to the saved baseline. The untethered enzyme is now labelled TadA8.20; numerical results are unchanged.
- Verified plotting in an isolated directory containing only the plotting script and its two CSV inputs, without sequencing summaries or the shared prepared table.
- Observed elapsed times, including container startup: preparation 6.76 seconds; plotting 7.49 seconds. Peak memory was approximately 155 MiB and 200 MiB, respectively. These are measurements from this run, not runtime guarantees.

## Changelog

- **2026-09-24:** Separated Figure 2B recorder preparation into Python and plotting into R; preserved numerical results and plot appearance apart from the requested TadA8.20 label; removed the bundled time course and its output, and cleaned up redundant filtered-recorder files locally. Other figure analyses and raw-read processing remain unchanged.
