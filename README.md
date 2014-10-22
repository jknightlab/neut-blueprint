---
title: Annotating Neutrophil eQTL with Blueprint Data
author: Peter Humburg
date: Tue 21 Oct 2014
---





# Introduction
The aim of this analysis is to annotate eQTL obtained from the analysis
of neutrophils with the results of bisulfite sequencing and ChIP-seq
for several histone modifications generated by the Blueprint consortium.
BED files with ChIP-seq peaks and hyper- and hypo-methylated regions
all neutrophil samples were obtained form the Blueprint ftp server.

All eQTL significant at an FDR threshold of 0.05 will be annotated with
the distance to the closest hyper/hypo methylated site as well as the distance
to the available histone marks. The subsequent analysis of these distances
focuses on eQTL of particular interest, e.g. the lead SNP for each gene or
all peak eQTL.



To be able to assess the significance of the proximity between eSNPs and
features of interest we also generate the corresponding background distributions
by considering the distance between the features and all tested SNPs.



# Methylation data
## Data processing
The BED files provided by Blueprint do not conform to the BED specification.
This causes difficulties when trying to process them with BED file parsers.
To remedy this the first five columns of each file are extracted prior
to further analysis.


The re-formatted BED files are imported into R for further processing with 
Bioconductor. 



## Defining regions of hyper- and hypo-methylation
Regions considered to be either hyper- or hypo-methylated are defined as follows.
A genomic region is considered to be hyper-methylated in neutrophils if it
meets the criteria for hyper-methylation used in the analysis of the Blueprint consortium
(average methylation of >0.75 and all CpGs in the regions have methylation >0.5)
in at least two neutrophil samplesand the region isn't hypo-methylated in any of 
the neutrophil samples. Hypo-methylated regions are defined accordingly.



![Length distribution of hyper- and hypo-methylated regions](figure/regionLength.png) 

Considering the length distribution of hyper- and hypo-methylated regions it is apparent
that a lot of regions are quite short while some, especially for hyper-methylation,
are quite large. This is further emphasised by the following summary table.


|    &nbsp;     |  hypo methylated  |  hyper methylated  |
|:-------------:|:-----------------:|:------------------:|
|   **Min.**    |         2         |         1          |
|  **1st Qu.**  |        54         |        894         |
|  **Median**   |        157        |        2715        |
|   **Mean**    |       455.6       |        4769        |
|  **3rd Qu.**  |        529        |        6367        |
|   **Max.**    |       24740       |       129600       |

Table: Summary of length distribution for hypo- and hyper-methylated regions.

Based on this we set the minimum overlap between regions from different individuals
to 50bp. 

```
## Warning: Removed 252 rows containing non-finite values (stat_density).
## Warning: Removed 45 rows containing non-finite values (stat_density).
```

![plot of chunk overlapRegions](figure/overlapRegions.png) 

While this excludes small regions the effect on the length of the remaining regions is limited,
consistent with substantial overlaps between samples. A total of 
39892 hypo- and
422555 hyper-methylated
regions remain after filtering. 


|    &nbsp;     |  hypo methylated  |  hyper methylated  |
|:-------------:|:-----------------:|:------------------:|
|   **Min.**    |        50         |         50         |
|  **1st Qu.**  |        250        |        1320        |
|  **Median**   |        618        |        3365        |
|   **Mean**    |       935.3       |        5470        |
|  **3rd Qu.**  |       1306        |        7297        |
|   **Max.**    |       22190       |       115400       |

Table: Summary of length distribution for hypo- and hyper-methylated consensus regions.

## Annotating eQTL
All significant eQTL are annotated with their distance to the closest hyper- and 
hypo-methylated regions.



Considering the most significant SNP for each gene, in addition to the
overal distribution of distances to the closest methylation island the
number of SNPs within islands, i.e. those with distance 0, are of
particular interest.

![Proportion of eSNPs contained within methylation islands](figure/methDistZero.png) 


The proportion of cis eSNPs located in hypo-methylated regions
(3.7365813%) is substantially
higher than would be expected by chance 
(only 0.8285134% 
of imputed SNPs are located in hypo-methylated regions). This corresponds
to an odd's ratio of 4.6. 
The significance of this difference is confirmed via Fisher's exact test.


|      P value      |  Alternative hypothesis  |
|:-----------------:|:------------------------:|
| _1.189e-59_ * * * |        two.sided         |

Table: Fisher's Exact Test for Count Data: 6933 and 181 out of 836800 and 4844 respectively.
  
For the much more common hyper-methylated sites the difference between eSNPs and
background is less pronounced. Overall 
42.7590822% of imputed SNPs
fall within a hyper-methylated region while only 
38.5631709% are located within one.
This corresponds to an odd's ratio of 
0.84. 
The significance of this difference is confirmed via Fisher's exact test.


|      P value      |  Alternative hypothesis  |
|:-----------------:|:------------------------:|
| _3.594e-09_ * * * |        two.sided         |

Table: Fisher's Exact Test for Count Data: 357808 and 1868 out of 836800 and 4844 respectively.
 
 The distribution of distances paint a similar picture. There is a pronunced
 shift in the location of eSNPs towards hypo-sensitive regions when compared
 to the distribution of all imputed SNPs. The difference observed in the distance
 from hyper-sensitive sites on the other hand is very small.
 
![Distance from lead SNP to nearest methylation island](figure/methDistLeadPlot.png) 

# Histone modification data
The Blueprint consortium provides BED files with peak calls for a number of 
histone modifications obtained from neutrophils. As with the methylation data
the format of the BED files is non-standard and the first five columns have to be
extracted prior to further processing.




There are 7 different histone modifications
in the data set with varying numbers of samples (see Table 1).


|     &nbsp;     |  Samples  |
|:--------------:|:---------:|
|  **H2A_Zac**   |     1     |
|  **H3K27ac**   |     7     |
|  **H3K27me3**  |     6     |
|  **H3K36me3**  |     7     |
|  **H3K4me1**   |     7     |
|  **H3K4me3**   |     7     |
|  **H3K9me3**   |     7     |

Table: **Table 1:** Summary of available ChIP-seq samples


# Appendix {-}
## Custom functions used

```r
loadBlueprint <- function(files) {
    splitNames <- strsplit(basename(files), ".", fixed = TRUE)
    sample <- sapply(splitNames, "[[", 1)
    type <- sapply(splitNames, "[[", 2)
    df <- data.frame(file = files, sample = sample, type = type)
    gr <- lapply(files, rtracklayer::import.bed, genome = "hg19", asRangedData = FALSE)
    list(meta = df, ranges = gr)
}

figRef <- local({
    tag <- numeric()
    created <- logical()
    used <- logical()
    function(label, caption, prefix = options("figcap.prefix"), sep = options("figcap.sep"), 
        prefix.highlight = options("figcap.prefix.highlight")) {
        i <- which(names(tag) == label)
        if (length(i) == 0) {
            i <- length(tag) + 1
            tag <<- c(tag, i)
            names(tag)[length(tag)] <<- label
            used <<- c(used, FALSE)
            names(used)[length(used)] <<- label
            created <<- c(created, FALSE)
            names(created)[length(created)] <<- label
        }
        if (!missing(caption)) {
            created[label] <<- TRUE
            paste0(prefix.highlight, prefix, " ", i, sep, prefix.highlight, 
                " ", caption)
        } else {
            used[label] <<- TRUE
            paste(prefix, tag[label])
        }
    }
})

tabRef <- local({
    tag <- numeric()
    created <- logical()
    used <- logical()
    function(label, caption, prefix = options("tabcap.prefix"), sep = options("tabcap.sep"), 
        prefix.highlight = options("tabcap.prefix.highlight")) {
        i <- which(names(tag) == label)
        if (length(i) == 0) {
            i <- length(tag) + 1
            tag <<- c(tag, i)
            names(tag)[length(tag)] <<- label
            used <<- c(used, FALSE)
            names(used)[length(used)] <<- label
            created <<- c(created, FALSE)
            names(created)[length(created)] <<- label
        }
        if (!missing(caption)) {
            created[label] <<- TRUE
            paste0(prefix.highlight, prefix, " ", i, sep, prefix.highlight, 
                " ", caption)
        } else {
            used[label] <<- TRUE
            paste(prefix, tag[label])
        }
    }
})
```

## Session Info

```
## R version 3.1.1 (2014-07-10)
## Platform: x86_64-pc-linux-gnu (64-bit)
## 
## locale:
##  [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8    
##  [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
##  [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] parallel  methods   stats     graphics  grDevices utils     datasets 
## [8] base     
## 
## other attached packages:
## [1] pander_0.3.8         scales_0.2.4         ggplot2_1.0.0       
## [4] GenomicRanges_1.16.4 GenomeInfoDb_1.0.2   IRanges_1.22.10     
## [7] BiocGenerics_0.10.0  knitr_1.6.14        
## 
## loaded via a namespace (and not attached):
##  [1] colorspace_1.2-4 digest_0.6.4     evaluate_0.5.5   formatR_1.0     
##  [5] grid_3.1.1       gtable_0.1.2     labeling_0.3     MASS_7.3-34     
##  [9] munsell_0.4.2    plyr_1.8.1       proto_0.3-10     Rcpp_0.11.2     
## [13] reshape2_1.4     stats4_3.1.1     stringr_0.6.2    tools_3.1.1     
## [17] XVector_0.4.0
```
