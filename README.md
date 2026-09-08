# PMOS trait operationalization: statistical analysis code

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22655208.svg)](https://doi.org/10.5281/zenodo.22655208)

Archived release DOI: [10.5281/zenodo.22655208](https://doi.org/10.5281/zenodo.22655208).
The concept DOI for all versions is [10.5281/zenodo.22655207](https://doi.org/10.5281/zenodo.22655207).

This repository contains the statistical analysis code for a cross-sectional
study of androgen, anti-Mullerian hormone (AMH), adiposity, and metabolic traits.
It covers NHANES August 2021-August 2023, NHANES 2017-March 2020, and a
local-only interface for the SWAN baseline public-use dataset.

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
4. evaluates alternative operational definitions and incremental fit;
5. runs weighted PCA, within-PSU parallel analysis, bootstrap stability, and
   permutation-based discordance benchmarks;
6. performs the earlier-NHANES comparison;
7. runs SWAN models only when a local file and variable mapping are available.

Generated files are written beneath `outputs/`, which is ignored by Git.

## SWAN local setup

Obtain ICPSR study 28762, version 5 under its applicable terms. Copy the
baseline data file to `data/restricted/` and edit
`config/swan_variable_map.csv` so each canonical variable points to the exact
column in the downloaded file. Set the environment variable `SWAN_FILE` to the
local file path before running. The script accepts `.dta`, `.sav`, `.sas7bdat`,
`.xpt`, `.rds`, or `.csv` files.

## Reproducibility notes

- Survey designs are created before domain restriction.
- Nonpositive or nonfinite subsample weights are excluded from sampled domains.
- Random procedures use fixed seeds from `config/analysis.yml`.
- The NHANES validation test checks the published cohort milestones of 824
  strictly eligible women and 643 women in the fully adjusted discovery model.
- SWAN code can be syntax-tested without the restricted file, but numerical
  reproduction requires the independently obtained dataset.

## Data access

NHANES data and documentation are available from the US National Center for
Health Statistics. SWAN baseline public-use data are available from ICPSR after
authentication and acceptance of the applicable terms. This repository does
not redistribute either dataset.

## License

Code is released under the MIT License. Dataset terms remain governed by NCHS
and ICPSR and are not altered by this software license.
