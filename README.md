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

Install the released version of mina from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("mina")
```
Or from github by:
```r
devtools::install_github("Guan06/MINA", dependencies = TRUE,
                          repos = c("https://cloud.r-project.org/",
                                    BiocManager::repositories()))
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
```r
hmp <- norm_tab(hmp, method = "total")

```
#### Composition-based diversity analysis
Dissimilarity / distance between pairwise samples were usually calculated by
comparing the differences indicated by the quantitative table and were used to
represent the beta-diversity of the community. Here we calculated weighted
Jaccard index
(https://en.wikipedia.org/wiki/Jaccard_index#Weighted_Jaccard_similarity_and_distance)
between pairwise HMP damples and Principal Coordinates Analysis (PCoA, also
referred as Classical multidimensional scaling MDS) was then used for
dimensionality reduction and visualization.
```r
hmp <- com_dis(hmp, method = "fJaccard")
# Unexplained variance ratio of the distance matrix, factors ordered according
# to the meta data shown in HMP website.
com_r2(hmp, group = c("sex", "RUNCENTER", "HMPBodysubsite", "area"))
hmp <- dmr(hmp)
p1 <- com_plot(hmp, match = "Sample_ID", color = "area")
p1
```
See full list of available distance by:
```r
?com_dis_list
```
Notably, we included TINA (Schmidt *et al.*, 2016) dissimilarity in the package.

#### Network-derived feature based diversit analysis
Afterwards, Spearman correlation between OTUs were calculated and coefficients
not less than 0.3 were retained in the sparse matrix for clustering by Markov
Cluster Algorithm (MCL, Enright, Dongen and Ouzounis, 2002) with parameter
'-I 2.5'. Later on, by summing up the abundance of OTUs belong to the same
cluster, network cluster quantitative table was obatained.
```r
hmp <- adj(hmp, method = "spearman")
hmp <- net_cls(hmp, method = "mcl", cutoff = 0.3)
hmp <- net_cls_tab(hmp)
# calculate community distance matrix based on network cluster table
hmp_nc <- hmp@hmp@cls_tab
hmp_nc_dis <- com_dis(hmp_nc, method = "fJaccard")
get_r2(hmp_nc_dis, hmp_des,
       group = c("sex", "RUNCENTER", "HMPBodysubsite", "area"))

hmp_dmr <- dmr(hmp_nc_dis)
p2 <- pcoa_plot(hmp_dmr, hmp_des, match = "Sample_ID", color = "area")
p2
```

### Community network comparison
We developed a bootstrap-permutation based method to test the significance of
network differences. By subsampling and bootstrap, true networks were
constructed from original dataset of each environment as shown below. By
randomly swapping the metadata of samples, permutated datasets were generated.
Networks of pseudo conditions were then inferred from permutated dataset.
Afterwards, true network distances (F) between each pairwise true networks and
pseudo distances (Fp) between each pairwise pseudo networks were calculated and
compared. An empirical P-value is then then calculated as
(Count Fp > F + 1) / (N + 1), where N is the total time of comparison between Fp
and F. Clearly, in order to observe a significant result, N need to be large
enough. By introducing both bootstrap and permutation process, the network
inference time is reduced to the sum of bootstrap and permutation time
(b1 + b2 + p1 + p2), resulting in retrench of computing time and space usage.

Here we use the maize data as example, networks of samples from different
compartments and host developmental stages were compared.
```r
maize <- norm_tab(maize, method = "raref", depth = 2000)
maize <- fit_tabs(maize)
maize <- bs_pm(maize, group = "Compartment", g_size = 200, s_size = 80)
maize <- net_dis(maize, method = "spectra")
maize@dis_stat

```
## References

