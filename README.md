# OrganAgeR

Organ Age Prediction from SomaScan Proteomic Data

## Installation
```r
# Install from GitHub
remotes::install_github("kaiyong821/OrganAgeRpackage")
```

## Usage
```r
library(OrganAgeR)

# Predict organ ages
result <- estimate_organ_age(
  data = your_data,
  age_col = "AGE",
  sex_col = "Sex_F",
  id_col = "HELPFul_ID"
)
```

## Requirements

- SomaScan 7K assay data
- Raw RFU values (NOT log-transformed)
- Exported from anmlSMP.ADat file
