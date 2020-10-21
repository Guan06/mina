# **M**icrobial community d**I**versity and **N**etwork **A**nalysis with **MINA**

An increasing number of microbiome datasets have been generated and analyzed with the help of rapidly developing sequencing technologies. At present, analysis of taxonomic profiling data is mainly conducted using composition-based methods, which ignores interactions between community members. Besides this, a lack of efficient ways to compare microbial interaction networks limited the study of community dynamics. To better understand how community diversity is affected by complex interactions between its members, we developed a framework (**M**icrobial community d**I**versity and **N**etwork **A**nalysis, **MINA**), a comprehensive framework for microbial community diversity analysis and network comparison.  By defining and integrating network-derived community features, we greatly reduce noise-to-signal ratio for diversity analyses.  A bootstrap and permutation-based method was implemented to assess community network dissimilarities and extract discriminative features in a statistically principled way.

## Overview of the workflow

The whole workflow could be divided into two main part: community diversity
analysis (green functions shown below) and network analysis (blue functions).
The `mina` object is shown here but one could also input pre-defined matrix
at each step (details in manual).

![whole_workflow](https://github.com/Guan06/MINA/blob/master/data-raw/workflow.png)

## Installation

The official version of `mina` can be installed from github by:
```r
devtools::install_github("Guan06/MINA", dependencies = TRUE,
                          repos = c("https://cloud.r-project.org/",
                                    BiocManager::repositories()))
```

## Overview of the workflow
The **MINA** workflow could be divided into two main parts: a) community diversity analysis (green functions shown below) and b) network analysis (blue functions). We define a data structure called mina object, which contains all relevant community features and can be used for every step in the analysis pipeline. Alternatively, the user can perform individual steps on pre-defined feature matrices (e.g. ASV / OTU tables) separately (see further details in the user manual).

**MINA** expects count data such as the commonly used OTU or ASV table to indicate the abundance of each community member in each sample. In addition, a descriptive metadata table is required for downstream analysis (e.g. comparison between treatments). Two example datasets were included in the package.  A detailed demonstration of the workflow, description of the data format, parameters and usage could be found in the accompanying vignette.

## Example of a diversity and network analysis with MINA

### Loading Input data
We included an OTU table of downloaded from [Human Microbiome Project](https://www.hmpdacc.org/hmp/HMQCP/) and an ASV table from the maize root microbiome (Bourceret and Guan *et al.*, *in prep*). These two datasets can be used to follow the entire **MINA** workflow with examples.
To import the data and create new `mina` object:
```r
hmp <- new("mina", tab = hmp_otu, des = hmp_des)
maize <- new("mina", tab = maize_asv, des = maize_des)
```
Altenatively, the two datasets are also included as `mina` objects in the package.

#### Data normalization
For community analyses, the data first needs to be normalized or rarefied before distance/dissimilarity matrix calculation.
```r
hmp <- norm_tab(hmp, method = "total")
```

### Composition-based diversity analysis
Pairwise dissimilarity/distance are usually calculated by comparing entries in a table of community features (e.g. species relative abundances). For instance, we can obtain a distance matrix for the HMP samples with the `com_dis` function:
```r
hmp <- com_dis(hmp, method = "bhjattacharyya")
```
The default method is for this calculation is Bray-Curtis, but the user can specify different distance/dissimilarity measures (see `?com_dis_list` for the full list). Once we have obtained pairwise distances, we can assess the amount of variance explained by our biological factors by running:
```r
com_r2(hmp, group = c("Sex", "Run_center", "Subsite", "Site"))
```
Next, we can perform a Principal Coordinates Analysis (PCoA, also referred as Classical multidimensional scaling MDS) for dimensionality reduction and visualization.
```r
hmp <- dmr(hmp)
com_plot(hmp, match = "Sample_ID", color = "Site")
```
<img src="https://github.com/Guan06/MINA/blob/master/data-raw/p1.png" alt="OTUs-based diversity" width="350" height="350">

### Diversity analysis based on network-derived features
To obtained network-derived features we first need to calculate an adjacency matrix based using the `adj` function after removing low prevalence community members to speed up the calculations.
```r
lst <- rownames(hmp@norm)[rowSums(hmp@norm > 0) > 50]
hmp@adj <- hmp@adj[lst, lst]dim(hmp@adj)
hmp <- adj(hmp, method = "spearman")
```
Next, we perform clustering on the community network using the Markov Cluster Algorithm ([MCL](https://micans.org/mcl/)) after removing weak edges in the adjacency matrix (<= 0.4 correlation coefficient) and obtain a sparse matrix representation of the community network. Both steps are performed simultaneously using the function `net_cls`.
```r
hmp <- net_cls(hmp, method = "mcl", cutoff = 0.4)
```
Once we have classified community members in terms of their interactions using the correlation network, we can obtain a network feature table by adding up the relative abundances of all members in each cluster:
```r
hmp <- net_cls_tab(hmp)
```
This new table represents the abundance of each network cluster in each sample. We can obtain network-based distance or dissimilarity matrices for diversity analyses as performed before:
```r
hmp_nc <- hmp@cls_tab
hmp_nc_dis <- com_dis(hmp_nc, method = "bhjattacharyya")
```
By calculating performing dimensionality reduction and calculating the percentage of unexplained variance using this feature table instead of the initial ASV/OTU table, we can see a marked increase in the signal-to-noise ratio:
```r
get_r2(hmp_nc_dis, hmp_des, group = c("Sex", "Run_center", "Subsite", "Site"))
hmp_dmr <- dmr(hmp_nc_dis)
pcoa_plot(hmp_dmr, hmp_des, match = "Sample_ID", color = "Site")
```
<img src = "https://github.com/Guan06/MINA/blob/master/data-raw/p2.png" alt = "Network Clusters-based diversity" width = 350 height = 350>

### Community network comparison
To test whether differences between community co-occurrence networks are statistically significant, we developed a bootstrap and permutation-based method.
By subsampling and bootstrap, true networks were constructed from original dataset of each environment as shown below. By randomly swapping the metadata of samples, permutated datasets were generated.
First, we use a subsampling approach to generate multiple bootstrap instances of each network, among which distances can be calculated. Next, we compare these distances with those inferred from bootstrap networks obtained after random permutation of sample labels.
An empirical P-value is then calculated by estimating how frequently the distance observed between true networks (*F*) is larger than the distance observed between permutated networks (*Fp*):
*P = (Count Fp > F + 1) / (N + 1)*
where *N* is the total time of comparison between *Fp* and *F*.
![workflow](https://github.com/Guan06/MINA/blob/master/data-raw/bs_pm.png)
Here we can use the dataset from the maize root microbiome as an example where we compare networks from samples obtained from different compartments and host developmental stages. As before, we first normalize the ASV table:
```r
maize <- norm_tab(maize, method = "raref", depth = 1000)
maize <- fit_tabs(maize)
```
Next, we can calculate the bootstrap networks by using `bs_pm` as follows:
```
maize <- bs_pm(maize, group = "Compartment", g_size = 100, s_size = 50)
```
Finally, we perform the statistical analysis with the function `net_dis`.
```r
maize <- net_dis(maize, method = "Jaccard")
```
We can then explore the results by accessing the `dis_stat` slot in the mina object:
```r
maize@dis_stat
```

| Compare                 | Distance_Mean | Distance_SD | Distance_PM_Mean | Distance_PM_SD | N    | p           |
|-------------------------|---------------|-------------|------------------|----------------|------|-------------|
| rhizosphere_rhizosphere | 0.6455   | 0.0285 | 0.6659      | 0.0968    | 1296 | 0.4271 |
| rhizosphere_root        | 0.9262   | 0.0017 | 0.7834      | 0.0458    | 1296 | 0.0008  |
| root_root               | 0.7112   | 0.0315 | 0.7306      | 0.0710    | 1296 | 0.5343 |
