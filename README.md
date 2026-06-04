
# ggebp: GGE Biplots with Mosaic Plots in R

This is a standalone R function (`ggebp.r`) developed to  construct GGE-biplots using `ggplot2`.

### The context
I created this function in 2015 to implement the GxE analysis techniques described by Laffont, Hanafi, and Wright (Crop Science 47:990-996; 53:2332-2341). I used it for MET data analysis at the time because no biplot tools using ggplot2 were available, at least to my knowledge.

![ggebp output screenshot](demo.png)

### Technical summary
This implementation was built from the ground up to understand the foundations of the GxE analysis. It performs Singular Value Decomposition (SVD) to partition the combined genotype and genotype-by-environment variance into principal components, and generates a customizable graphic using `ggplot()` function. It also includes an interesting feature to generate a mosaic plot of the TSS partition (this was basically "stolen" from Laffont et al., 2007).

Legacy code from 2015, provided "as-is". It was written in base R, old-school style, prior to the widespread use of pipes and before the Tidyverse became the standard for data manipulation. I found it on an old drive and made it available to anyone interested. It’s still fully functional in 2026. 

### Usage
The function is designed to be user-friendly, even offering a file selection window if no data is provided.
```r
source("ggebp.r")

# Yang barley dataset included in the repo, also available in 'agridat' package
yb <- read.csv("yang_barley.csv")

# Basic plot
ggebp(yb, title = "Example basic Biplot")

# Change some plot parameters, includes SS partition (mosaic plot)
ggebp(yb, mosaic = T, obs.labels = T,
    var.factor = 1.5, obs.color="green3", var.color="multi", line.width = 1, 
    title = "Example improved Biplot + SSP Mosaic Plot")

```

