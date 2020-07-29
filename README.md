# Microbial community dIversity and Network Analysis with MINA

An increasing number of microbiome datasets have been generated and analysed
with the rapidly developing sequencing technologies. At present, analysis of
taxonomic profiling data is mainly conducted using composition-based methods,
which ignores interactions between community members. Besides, the
comprehensive comparative network analysis is unexplored, limiting the study
of community dynamics. The goal of *MINA* is to provide a thorough framework for
microbial community analysis based on higher order community features to
better understand the principles that govern the establishment of those
communities. We reduced the noise / signal ratio for diversity analysis by
integrating the network-derived features and introduced a
bootstrap-permutation based network comparison method to statistically assess
community networks dissimilarities under specific condition and to extract
discriminative features.

## Overview of the workflow

## Installation

You can install the released version of mina from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("mina")
```
## Introduction
Microbes play an important role in most of ecosystems and by interacting with
each other, they assemble into complex systems. The formation and stability of
microbial community is affected by both biotic and abiotic environmental factors.
The microbial components could be determined by marker gene and the development
of sequencing technology makes the high-throughput profiling of microbial
communities possible. Typically, analysis of this data includes estimating
within and between sample diversities (alpha- and beta-diversity, respectively)
based on compositions extracted from sequences, such as OTUs or ASVs, which are
operational taxonomic units clustered using arbitrary threshold or exact
sequence variants could be identified by single nucleotide difference.
By counting the number of observed compositions, distance or dissimilarity
between samples calculated from counts differentiation of compositions is used
to indicate the beta diversity of community.

Although compositional approaches provide a way to characterise community
structure and to measure the differences between samples, they capture only
static fetures and ignore the dynamics of the system, such as interactions
between compositions of the community. To overcome these limitations,
co-occurrence networks are typically inferred. In these microbial community
networks, nodes represent the community members and relationships between
microbes are indicated by undirected edges, which are inferred by comparing
the covariance of microbes across samples. The number of samples, therefore,
has an impact on the robustness of the constructed network. Typically, for
network comparison, the distance between correlation matrices are calculated
and defined as differences between networks, whereas the corresponding
statistical test methods are missing due to the computational limitation.

To better understand the assembly and maintenance of the community structure,
we developed a framework (Microbial community dIversity and Network Analysis,
MINA) for microbial community data processing. We implemented both
composition and netowrk derived feature based diversity analysis and in
addition, a bootstrap-permutation method is introduced to thoroughly compare
ecological networks and to assess their dissimilarity.

## Application and results
*MINA* expects count data such as common used OTU or ASV table to indicate the
abundance of each composition in each sample. Besides, a descriptive table
which indicates the meta data is required for later comparison analysis. Two
datasets were included in the package as demonstration and the details about
data format, parameters and usage could be found in vignette.

### Input data
We included OTU table of Human Microbiome Project downloaded from
https://www.hmpdacc.org/hmp/HMQCP/ and maize root-associated microbial
community profiling ASV table from Bourceret and Guan *et al*., 2020.
To import the data and create new 'mina' object:
```r
hmp <- new("mina", tab = hmp_otu, des = hmp_des)
maize <- new("mina", tab = maize_asv, des = maize_des)
```
### Community diversity analysis
For the community diversity analysis, here we use hmp data as example. The data
needed to be normalized / rarefied before distance / dissimilarity matrix
calculation.

### Community network comparison
