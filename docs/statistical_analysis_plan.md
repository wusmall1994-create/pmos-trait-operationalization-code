# Statistical analysis plan implemented by this release

## Study populations

- Discovery: NHANES August 2021-August 2023, women aged 20-44 years.
- Temporal comparison: NHANES 2017-March 2020 pre-pandemic, with the same age
  and reproductive eligibility rules.
- Cross-cohort comparison: SWAN baseline public-use data obtained separately
  from ICPSR.

NHANES exclusions are pregnancy or indeterminate pregnancy status, current
breastfeeding, hysterectomy, and bilateral oophorectomy. Designs are created
before domain restriction, and nonpositive subsample weights are excluded.

## Trait construction

Hormones are natural-log transformed, residualized for age using weighted
linear regression within period, and standardized with the relevant survey
weights. The androgen score is the mean of standardized total testosterone,
androstenedione, and DHEAS and requires at least two observed components. AMH is
analyzed separately. FAI equals 100 multiplied by testosterone in nmol/L and
divided by SHBG in nmol/L.

The nonfasting metabolic score is the mean of standardized HbA1c, inverse HDL-C,
and mean arterial pressure. The fasting score additionally includes standardized
log HOMA-IR and log triglycerides. HOMA-IR is fasting glucose in mg/dL multiplied
by fasting insulin in micro-international units/mL and divided by 405. Mean
arterial pressure is (systolic + 2 x diastolic) / 3.

## Models

Primary survey-weighted linear models use a fixed complete-case sample across
the adjustment sequence: exposure only; age and race/ethnicity; poverty-income
ratio and current smoking; waist-to-height ratio; and an
exposure-by-waist-to-height-ratio interaction. Estimates are expressed in
outcome SD per exposure SD.

The exploratory common-sample matrix crosses 10 hormone representations with
9 metabolic outcomes, before and after adiposity adjustment. Benjamini-Hochberg
control is applied within each matrix. Nested weighted R-squared analyses enter
three androgens, AMH, and inverse SHBG in both orders.

Weighted PCA uses a survey-weighted correlation matrix. Parallel analysis uses
500 independent within-PSU permutations per variable. Sampling stability uses
1,000 PSU-within-stratum bootstrap replicates with axis matching and sign
alignment. Directional discordance is compared with analytical independence and
1,000 within-PSU permutations.

Sensitivity analyses exclude all recent prescription users, restrict age to
20-39 years, exclude diabetes, remove extreme hormone values, substitute BMI or
waist circumference for waist-to-height ratio, evaluate FAI and SHBG-related
constructs, and omit one PSU at a time.

## Randomness and multiplicity

Seeds and replication counts are fixed in `config/analysis.yml`. Holm correction
is used for prespecified primary tests, and Benjamini-Hochberg correction is used
for exploratory families.

