
# ggebp: GGE Biplots with Mosaic Plots in R

This repository contains a standalone R function (`ggebp.R`) developed to perform Singular Value Decomposition (SVD) and construct GGE-biplots using `ggplot2`.

### The Context
I developed this tool back in 2015 to implement the advanced GxE analysis techniques described by Laffont, Hanafi, and Wright. The function integrates numerical and graphical measures from their 2007 work (Crop Science 47:990-996) and the specific GGEBP methodology from their 2013 paper (Crop Science 53:2332-2341).

I’m sharing this for anyone interested in implementing statistical theory into functional and practical code.

### Technical Highlights
Unlike standard "black-box" packages, this implementation was built from the ground up to provide transparency in the GxE (Genotype by Environment) analysis:

*   **Custom SVD Core:** Uses the `svd()` function for precise control over axis scaling and rotation.
*   **Variance Partitioning:** It explicitly calculates and prints the Sum of Squares (SS) partition for Genotype (SSG) and Interaction (SSGE).
*   **Mosaic Plot:** Includes a feature to generate a mosaic plot of the Total Sum of Squares (TSS) partition, showing the anatomy of variance at a glance.
*   **Automation:** Features an "auto-flip" logic that ensures the genotype ordinate is positively correlated with genotype means, saving the researcher from manual axis rotation.

### Usage
The function is designed to be user-friendly, even offering a file selection window if no data is provided.
```r
source("ggebp.R")

# Example:
# ggebp(my_data, mat = TRUE, mosaic = TRUE, title = "Yield Trials")

```

### Legacy Status

This code is provided "as-is" from my 2015 archives. It was written in Base R, before the Tidyverse became the standard for data manipulation. It reflects my early approach to solving complex problems using R's core functional programming capabilities.
