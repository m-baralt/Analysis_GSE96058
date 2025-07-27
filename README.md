# Analysis GSE96058

This repository contains code used to analyze the publicly available breast cancer dataset **GSE96058**.

## Contents

- `analysis_functions.R`: Contains R functions used in the different steps of the analysis pipeline.
- `Data_handling.R`: Downloads and prepares the dataset, runs the pipeline using the defined functions, and saves the results.
- `slides.qmd`: A Quarto presentation that loads the saved results and generates slides summarizing the analysis.

## Usage

1. Open `Data_handling.R` and run the code to process and analyze the data.
2. Once the results are saved, render `slides.qmd` to generate the presentation.

## Requirements
- R (≥ 4.0)
- Packages: `GEOquery`, `magrittr`, `ggplot2`, `RSpectra`, `plyr`, `survival`, `survminer`, `gprofiler2`, `doParallel`, `ranger`
