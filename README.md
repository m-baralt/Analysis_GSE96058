# Analysis GSE96058

This repository contains code used to analyze the publicly available breast cancer dataset **GSE96058**.

## Contents

- `analysis_functions.R`: Contains R functions used in the different steps of the analysis pipeline.
- `Data_handling.R`: Downloads and prepares the dataset, runs the pipeline using the defined functions, and saves the results.
- `slides.qmd`: A Quarto presentation that loads the saved results and generates slides summarizing the analysis.

## Usage

1. **Clone this repository:**

```
git clone https://github.com/m-baralt/Analysis_GSE96058.git
```

2. Navigate to the repository directory.

3. Open the project in RStudio by double-clicking the file: `Analysis_GSE96058.Rproj`. This ensures the working directory is correctly set, so all relative paths will work automatically.

4. Download the expression dataset using the following command lines:

```
chmod +x download_expression_data.sh
./download_expression_data.sh
```

This will download the expression dataset into the `data/` folder.

5. Open `Data_handling.R` and run the code to process and analyze the data.

6. Render the `slides.qmd` file to generate the presentation after results are saved.

*OR*

If you prefer, you can skip step 5 by downloading precomputed results using:

```
chmod +x download_classical_statistics.sh
chmod +x download_survivalRF.sh
./download_classical_statistics.sh
./download_survivalRF.sh
```

Then render the `slides.qmd` file from Rstudio.

## Requirements
- R (≥ 4.0)
- Packages: `GEOquery`, `magrittr`, `ggplot2`, `RSpectra`, `plyr`, `survival`, `survminer`, `gprofiler2`, `doParallel`, `ranger`
