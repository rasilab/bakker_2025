# BoxB Figure Generation Script Modularization Plan

## Executive Summary

The `analysis/boxb_in_vitro_sequencing/scripts/generate_figures.R` script is a large monolithic file (1,129 lines) that performs data loading, statistical analysis, and figure generation for BoxB hairpin *in vitro* sequencing analysis. This plan outlines how to split it based on **figure-specific data dependencies** into 5 focused modules that can run independently and quickly.

## Current Script Analysis

**Total Lines:** 1,129
**Purpose:** Analyze *in vitro* sequencing data for BoxB hairpins, calculate editing fractions, and generate publication figures
**Dependencies:** tidyverse, plyranges, rasilabRtemplates, ggpubr, RColorBrewer, Biostrings, cowplot

## Figure-to-Data Dependency Analysis

### Core Data Sources
- **`summary_stats_combined/*.csv.gz`** - Primary editing statistics for recorder sequences
- **`boxb_stats_combined/*.csv.gz`** - BoxB loop-specific editing statistics  
- **`barcode_annotations.csv`** + **`sample_info.csv`** - Required by all analyses
- **`hairpin_annotations.tsv`** - BoxB variant sequences and predicted stability

### Figure Groupings by Data Dependencies

**Group 1: Recorder-Only Analysis**
- Uses: `summary_stats_combined` + core annotations
- Figures: `figure_2b` (recorder part), `boxb_distance.pdf`, `figure_2d_2e.pdf`, context plots

**Group 2: Loop-Only Analysis** 
- Uses: `boxb_stats_combined` + core annotations
- Figures: `figure_2b` (loop part)

**Group 3: BoxB Variant Analysis**
- Uses: `summary_stats_combined` + core annotations + `hairpin_annotations.tsv`
- Figures: `gnra.pdf`, `figure_3c_3d.pdf`, `figure_3i_3j.pdf`, supplementary figures

**Group 4: Multi-Source Comparisons**
- Uses: Both `summary_stats_combined` + `boxb_stats_combined` + all annotations
- Figures: `total_edits_types.pdf`, recruitment comparisons, correlation plots

## Proposed Modular Structure

### Module 1: Data Preparation
**File:** `01_prepare_data.R`
**Complexity:** Low-Medium
**Runtime:** 2-5 minutes (I/O bound)

**Purpose:**
- Load and merge all raw data files once
- Create standardized, analysis-ready data frames
- Save processed data as `.rds` files for fast loading

**Inputs:**
- `../data/summary_stats_combined/*.csv.gz`
- `../data/boxb_stats_combined/*.csv.gz`  
- `../annotations/barcode_annotations.csv`
- `../annotations/sample_info.csv`
- `../annotations/hairpin_annotations.tsv`

**Outputs:**
- `../tables/target_data.csv` (recorder editing data with annotations)
- `../tables/loop_data.csv` (BoxB loop editing data with annotations)
- `../tables/hairpin_annotations.csv` (processed hairpin metadata)

### Module 2: Figure 2 Generation (Concentration & Distance Effects)
**File:** `02_generate_figure2.R`
**Complexity:** Medium
**Runtime:** 30-60 seconds

**Purpose:**
- Generate Figure 2 panels showing concentration-response and distance effects
- Focus on core editing efficiency measurements

**Inputs:**
- `../tables/target_data.csv`
- `../tables/loop_data.csv`

**Outputs:**
- `../figures/figure_2b.pdf` (concentration-response curves)
- `../figures/boxb_distance.pdf` (distance dependence)
- `../figures/figure_2d_2e.pdf` (sequence context effects)
- `../tables/mean_editing_per_concentration.tsv`
- `../tables/mean_editing_recorder_position.tsv`

**Key Analysis:** Wild-type controls, concentration-response curves, distance-decay analysis

### Module 3: Figure 3 Generation (BoxB Variants)
**File:** `03_generate_figure3.R`
**Complexity:** Medium
**Runtime:** 30-60 seconds

**Purpose:**
- Generate Figure 3 panels analyzing BoxB loop and stem variants
- Focus on sequence-structure relationships

**Inputs:**
- `../tables/target_data.csv`
- `../tables/hairpin_annotations.csv`

**Outputs:**
- `../figures/gnra.pdf` (GNRA motif analysis)
- `../figures/figure_3c_3d.pdf` (loop variant effects)
- `../figures/figure_3i_3j.pdf` (stem stability effects)
- `../tables/editing_per_loop_variant.tsv`
- `../tables/mean_editing_per_stem_variant.tsv`

**Key Analysis:** GNRA motif effects, base-pairing analysis, free energy correlations

### Module 4: Figure 4 & Recruitment Analysis
**File:** `04_generate_figure4_recruitment.R`
**Complexity:** Medium-High
**Runtime:** 45-90 seconds

**Purpose:**
- Generate Figure 4 and recruitment strategy comparisons
- Cross-platform analysis requiring multiple data sources

**Inputs:**
- `../tables/target_data.csv`
- `../tables/loop_data.csv`
- `../tables/hairpin_annotations.csv`

**Outputs:**
- `../figures/total_edits_types.pdf`
- `../figures/gnra_recruitment_strategies.pdf`
- `../figures/recruitment_heatmaps_8_10.pdf`
- `../figures/recruitment_stem_stability.pdf`
- `../figures/tada_vs_ln_tada.pdf` (correlation plots)
- `../tables/mean_editing_per_recruitment_type.tsv`

**Key Analysis:** Multi-platform comparisons, recruitment strategy effects, correlation analysis

### Module 5: Supplementary Figures
**File:** `05_generate_supplementary_figures.R`
**Complexity:** Medium
**Runtime:** 30-60 seconds

**Purpose:**
- Generate all supplementary figures and supporting analyses
- Extended context analysis and validation plots

**Inputs:**
- `../tables/target_data.csv`
- `../tables/hairpin_annotations.csv`

**Outputs:**
- `../figures/sup_fig_1a_1b.pdf`
- `../figures/figure_S2a.pdf`, `../figures/figure_S2b.pdf`
- `../figures/figS2_c_d.pdf`
- `../figures/figure_S3_c.pdf`
- `../figures/context_constant_recruitments.pdf`
- Additional context analysis tables

**Key Analysis:** Extended sequence context, validation controls, supplementary comparisons

## Implementation Strategy

### Phase 1: Extract Data Preparation (Module 1)
1. Create `01_prepare_data.R` with lines 1-75 from original script
2. Add `write_csv()` calls for `target_data`, `loop_data`, `hairpin_annotations`
3. Test data loading independence and validate merged outputs
4. Ensure `../tables/` directory exists for output files

### Phase 2: Extract Figure-Specific Modules (Modules 2-5)
1. **Module 2:** Extract Figure 2 generation code, replace data loading with `read_csv()`
2. **Module 3:** Extract Figure 3 generation code, add hairpin annotation processing
3. **Module 4:** Extract Figure 4 + recruitment analysis, handle multi-source dependencies
4. **Module 5:** Extract supplementary figure code, ensure independence from main figures

### Phase 3: Validation & Integration
1. Create `run_all_modules.R` master script
2. Validate that all outputs match original script exactly
3. Add error handling and progress reporting
4. Test module independence (ability to run individual modules)

## Dependencies and Execution Patterns

```
Module 1 (Data Preparation)
    ↓
Module 2, 3, 4, 5 (Figure Generation) - ALL can run in parallel
```

**Independent Execution Examples:**
- Update Figure 2 aesthetics: Run only Module 2 (30 seconds)
- Analyze new BoxB variants: Run Module 1 + Module 3 (3-4 minutes)
- Generate recruitment comparison: Run Module 1 + Module 4 (3-5 minutes)

## Benefits of Figure-Based Modularization

1. **Targeted Development:** Modify specific figure aesthetics without touching other analyses
2. **Data Efficiency:** Each module loads only the data it needs
3. **Debugging:** Isolate issues to specific figure generation logic
4. **Parallel Execution:** All figure modules can run simultaneously after data prep
5. **Incremental Analysis:** Add new figures without affecting existing ones
6. **Publication Workflow:** Generate specific figures for manuscript revisions quickly

## Estimated Runtimes (After Modularization)

- **Module 1 (Data Prep):** 2-5 minutes (one-time cost)
- **Module 2 (Figure 2):** 30-60 seconds
- **Module 3 (Figure 3):** 30-60 seconds  
- **Module 4 (Figure 4):** 45-90 seconds
- **Module 5 (Supplementary):** 30-60 seconds

**Total runtime:** ~5-10 minutes (when running all modules)
**Figure update runtime:** 30-90 seconds (when updating specific figures)
**Parallel runtime:** ~3-7 minutes (when running figure modules in parallel)

## Implementation Notes

- Preserve exact statistical calculations and figure outputs
- Use consistent theming and color schemes across modules  
- Include module-specific error handling and progress indicators
- Add command-line arguments for figure selection
- Document data dependencies clearly in each module
- Test independence by running modules out of order