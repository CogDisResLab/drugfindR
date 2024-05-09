## v0.99.596 (2024-05-09)

## v0.99.595 (2024-05-09)

### Feat

- have the vignette load the data from the local file
- **vignette**: change output to Bioconductor format
- **vignette**: add the data file locally

### Refactor

- moved to camelCase naming

## v0.99.583 (2024-04-22)

### Feat

- fix the filtering criteria for signature

### Fix

- **get_signature**: remove the redundant "Error" in the message of `stop()`
- fix a building error for investigate_signature
- **get_signature**: remove the redundant "Error" in the message of `stop()`

## v0.99.559 (2023-09-19)

## v0.99.558 (2023-09-14)

## v0.99.557 (2023-09-14)

## v0.99.556 (2023-09-14)

### BREAKING CHANGE

- closes #14

### Feat

- **get_concordants**: return error on API connection failures

### Fix

- return a more informative error on `get_signature`

## v0.99.552 (2023-09-14)

### BREAKING CHANGE

- closes #21

### Feat

- **get_signature**: get_signature returns an error instead of silently proceeding

## v0.99.551 (2023-09-14)

### Feat

- update the parameter name
- bundle commit for linting and styling
- pass preliminary checks

### Fix

- **prepare_signature**: fix the warnings for using the .data pronoun

### Refactor

- **raw-data**: compress raw data to save space
