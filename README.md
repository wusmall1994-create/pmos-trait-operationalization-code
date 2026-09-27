# PMOS trait operationalization: statistical analysis code

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22847160.svg)](https://doi.org/10.5281/zenodo.22847160)

Archived release DOI: [10.5281/zenodo.22847160](https://doi.org/10.5281/zenodo.22847160).
The concept DOI for all versions is [10.5281/zenodo.22655207](https://doi.org/10.5281/zenodo.22655207).

This repository contains the statistical analysis code for cross-sectional and
longitudinal analyses of sex hormone-binding globulin (SHBG), androgen,
anti-Mullerian hormone (AMH), adiposity, and metabolic traits. It covers NHANES
August 2021-August 2023, NHANES 2017-March 2020, and repeated SWAN visits.

## Code-only release

The repository intentionally excludes:

- manuscript and submission files;
- NHANES or SWAN participant-level data;
- generated tables, figures, and model outputs;
- access tokens, local paths, email addresses, telephone numbers, and other
  personal contact information.

NHANES files are downloaded at run time from official NCHS endpoints. The SWAN
public-use package must be obtained independently from ICPSR and placed locally;
the package is never copied into version control.

## Requirements

- R 4.3 or later (the reported analysis used R 4.6.1)
- PowerShell 7 or Windows PowerShell for `run_all.ps1` (optional)
- R packages listed in `DESCRIPTION`

Install dependencies with:

```r
source("R/00_dependencies.R")
install_dependencies()
```

## Run

From the repository root:

```powershell
./run_all.ps1
```

or:

```r
source("run_all.R")
```

The pipeline performs the following steps:

1. downloads and checksum-records official NHANES XPT files;
2. constructs the eligible analytic domains and derived traits;
3. fits prespecified survey-weighted models and sensitivity analyses;
4. evaluates alternative operational definitions, replays eight
   literature-aligned biochemical androgen definitions, and estimates
   incremental fit;
5. runs weighted PCA, within-PSU parallel analysis, bootstrap stability, and
   permutation-based discordance benchmarks;
6. performs the earlier-NHANES comparison;
7. runs the optional SWAN baseline analysis when `SWAN_FILE` and its variable
   mapping are available;
8. runs optional longitudinal SWAN mixed models, mutually adjusted models,
   medication and menopause-stage sensitivity analyses, continuous-time AR(1)
   models, and next-visit models when `SWAN_DIR` is available;
9. generates the aggregate analyses for Supplementary Tables S32-S37,
   including extended adjustment, cohort comparison, missingness, participant
   flow, visit contributions, and the two-fold SHBG interpretation.

Generated files are written beneath `outputs/`, which is ignored by Git.

## SWAN baseline setup

Obtain ICPSR study 28762, version 5 under its applicable terms. Copy the
baseline data file to `data/restricted/` and edit
`config/swan_variable_map.csv` so each canonical variable points to the exact
column in the downloaded file. Set the environment variable `SWAN_FILE` to the
local file path before running. The script accepts `.dta`, `.sav`, `.sas7bdat`,
`.xpt`, `.rds`, or `.csv` files.

## SWAN longitudinal setup

Obtain the required SWAN public-use visit files directly from ICPSR under the
applicable terms. Keep all `*-Data.dta` files outside version control and set
`SWAN_DIR` to the directory containing them. From the repository root, run:

```powershell
Rscript --vanilla analysis/assess_swan_longitudinal.R
Rscript --vanilla analysis/run_swan_longitudinal_analysis.R
```

The assessment script reports visit-level completeness and repeat-measurement
eligibility. The longitudinal script performs visit-specific standardization,
within-between decomposition, mixed-effects analyses, adiposity adjustment,
mutually adjusted hormone models, sensitivity analyses, and exploratory
next-visit models. It then invokes `analysis/run_swan_reporting_extensions.R`
to generate Supplementary Tables S32-S37. Only aggregate tables and figures
are written beneath the ignored `outputs/` directory; participant-level
analytic data are not exported.

The study-design figure can be regenerated with
`analysis/build_study_design_figure.R`. Set `PMOS_FIGURE_DIR` to override its
default destination under `outputs/figures`.

## Reproducibility notes

- Survey designs are created before domain restriction.
- Nonpositive or nonfinite subsample weights are excluded from sampled domains.
- Random procedures use fixed seeds from `config/analysis.yml`.
- The NHANES validation test checks the published cohort milestones of 824
  strictly eligible women and 643 women in the fully adjusted discovery model.
- SWAN code can be syntax-tested without restricted files, but numerical
  reproduction requires the independently obtained visit datasets.

## Data access

NHANES data and documentation are available from the US National Center for
Health Statistics. SWAN baseline public-use data are available from ICPSR after
authentication and acceptance of the applicable terms. This repository does
not redistribute either dataset.

## License

Code is released under the MIT License. Dataset terms remain governed by NCHS
and ICPSR and are not altered by this software license.
