# **M**icrobial community d**I**versity and **N**etwork **A**nalysis with **MINA**

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

The whole workflow could be divided into two main part: community diversity
analysis (green functions shown below) and network analysis (blue functions).
The `mina` object is shown here but one could also input pre-defined matrix
at each step (details in manual).

![whole_workflow](https://github.com/Guan06/MINA/blob/master/data-raw/workflow.png)

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
we developed a framework (**M**icrobial community d**I**versity and **N**etwork **A**nalysis,
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
We included OTU table of downloaded from [Human Microbiome Project](https://www.hmpdacc.org/hmp/HMQCP/)
and maize root-associated microbial community profiling ASV table (unpublished)
In the HMP dataset, 2711 samples with 27,627 OTUs and for maize dataset, 437
samples with 11,098 ASVs were presented.
To import the data and create new `mina` object:
```r
hmp <- new("mina", tab = hmp_otu, des = hmp_des)
maize <- new("mina", tab = maize_asv, des = maize_des)
```
This step could also be skipped since the two dataset are included as 'mina'
object in the package.
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
represent the beta-diversity of the community. Here we calculated `bhjattacharyya`
between pairwise HMP damples and Principal Coordinates Analysis (PCoA, also
referred as Classical multidimensional scaling MDS) was then used for
dimensionality reduction and visualization.
```r
# Distance matrix calculation, took around 10 min for a computer cluster with 48
# threads and 756 Gb memory.
hmp <- com_dis(hmp, method = "bhjattacharyya")
# Unexplained variance ratio of the distance matrix, factors ordered according
# to the meta data shown in HMP website.
# R2 = 0.644
com_r2(hmp, group = c("Sex", "Run_center", "Subsite", "Site"))

hmp <- dmr(hmp)
p1 <- com_plot(hmp, match = "Sample_ID", color = "Site")
```
The OTUs-based community diversities:
```r
p1
```

<img src="https://github.com/Guan06/MINA/blob/master/data-raw/p1.png" alt="OTUs-based diversity" width="350" height="350">

See full list of available distance by:
```r
?com_dis_list
```
Notably, we included TINA ([Schmidt *et al.*, 2016](https://doi.org/10.1038/ismej.2016.139)) dissimilarity in the package.

#### Network-derived feature based diversit analysis
Afterwards, Spearman correlation between OTUs were calculated and coefficients
not less than 0.4 were retained in the sparse matrix for clustering by Markov
Cluster Algorithm ([MCL](https://micans.org/mcl/)) with parameter
'-I 2.5'. Later on, by summing up the abundance of OTUs belong to the same
cluster, network cluster quantitative table was obatained. Since there are much
less clusters than compositions, the computing time would be decreased
dramatically as well.
```r
# Adjacency matrix calculation, 12632 seconds used.
hmp <- adj(hmp, method = "spearman")

# Remove OTUs appeared in not more than 50 samples, 10577 OTUs were remained
lst <- rownames(hmp@norm)[rowSums(hmp@norm > 0) > 50]
hmp@adj <- hmp@adj[lst, lst]
dim(hmp@adj)

# 2166 components are removed for clustering because no strong edge was found
# between those OTUs.
hmp <- net_cls(hmp, method = "mcl", cutoff = 0.4)
# 13069 seconds used for this step, ~3.63 hours

# 800 clusters containing 8411 OUTs
hmp <- net_cls_tab(hmp)

# calculate community distance matrix based on network cluster table
hmp_nc <- hmp@cls_tab
hmp_nc_dis <- com_dis(hmp_nc, method = "bhjattacharyya")
get_r2(hmp_nc_dis, hmp_des,
       group = c("Sex", "Run_center", "Subsite", "Site"))
# 0.399

hmp_dmr <- dmr(hmp_nc_dis)
p2 <- pcoa_plot(hmp_dmr, hmp_des, match = "Sample_ID", color = "Site")
```
So by integrating the abundances of closely related OTUs to the abundances of
clusters, the unexplained variance ratio decrease from 0.644 to 0.399.

The network clusters-based community diversities:

<img src = "https://github.com/Guan06/MINA/blob/master/data-raw/p2.png" alt = "Network Clusters-based diversity" width = 350 height = 350>

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

![workflow](https://github.com/Guan06/MINA/blob/master/data-raw/bs_pm.png)

Here we use the maize data as example, networks of samples from different
compartments and host developmental stages were compared. After filtering,
rhizosphere and root networks were infered and compared.

```r
maize <- norm_tab(maize, method = "raref", depth = 1000)
maize <- fit_tabs(maize)

# 10851 components are used for bs_pm before filtering.
# 277 samples are used for bs_pm before filtering.
# 2 groups with 237 samples used for bootstrap.
# 192 seconds were used for this step
maize <- bs_pm(maize, group = "Compartment", g_size = 100, s_size = 50)
maize <- net_dis(maize, method = "Jaccard")
```
And the distance statistic results:
```r
maize@dis_stat
```

| Compare                 | Distance_Mean | Distance_SD | Distance_PM_Mean | Distance_PM_SD | N    | p           |
|-------------------------|---------------|-------------|------------------|----------------|------|-------------|
| rhizosphere_rhizosphere | 0.6455   | 0.0285 | 0.6659      | 0.0968    | 1296 | 0.4271 |
| rhizosphere_root        | 0.9262   | 0.0017 | 0.7834      | 0.0458    | 1296 | 0.0008  |
| root_root               | 0.7112   | 0.0315 | 0.7306      | 0.0710    | 1296 | 0.5343 |

Since the bootstrap and permutation, the distance calculation process is
nondeterministic, however the result and conclusion should not change a lot.

## Conclusion
The package provides thorough microbial diversity and network analysis tools.
Network cluster based diversity analysis reduced the noise.
By introducing the bs-pm process, we decreased the demonding time, space and computing
capacity dramatically.
