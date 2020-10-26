pkgname <- "mina"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('mina')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("adj-matrix")
### * adj-matrix

flush(stderr()); flush(stdout())

### Name: adj,matrix,ANY-method
### Title: Calculate the adjacency matrix of @norm by correlation with
###   matrix as input.
### Aliases: adj,matrix,ANY-method adj,matrix,character-method

### ** Examples

asv_norm <- norm_tab(maize_asv2, method = "raref", depth = 1000)
asv_adj <- adj(asv_norm, method = "pearson")



cleanEx()
nameEx("adj-mina")
### * adj-mina

flush(stderr()); flush(stdout())

### Name: adj,mina,ANY-method
### Title: Calculate the adjacency matrix of @norm by correlation with
###   'mina' class object as input.
### Aliases: adj,mina,ANY-method adj,mina,character-method

### ** Examples

maize <- new("mina", tab = maize_asv2, des = maize_des2)
maize <- norm_tab(maize, method = "raref")
maize <- fit_tabs(maize)
maize <- adj(maize, method = "spearman", sig = FALSE)



cleanEx()
nameEx("adj")
### * adj

flush(stderr()); flush(stdout())

### Name: adj
### Title: Calculate the adjacacency matrix of @norm by correlation.
### Aliases: adj

### ** Examples

maize <- new("mina", tab = maize_asv2, des = maize_des2)
maize <- norm_tab(maize, method = "raref", depth = 1000)
maize <- fit_tabs(maize)
maize <- adj(maize, method = "spearman", sig = FALSE)



cleanEx()
nameEx("adj_method_list")
### * adj_method_list

flush(stderr()); flush(stdout())

### Name: adj_method_list
### Title: List of adjacency matix calculation methods/ orrelations
###   supported in 'adj'
### Aliases: adj_method_list
### Keywords: datasets

### ** Examples

? adj_method_list



cleanEx()
nameEx("bs_pm-mina")
### * bs_pm-mina

flush(stderr()); flush(stdout())

### Name: bs_pm,mina,ANY-method
### Title: Inferring the network of different group of samples and test
###   significance by permutation.
### Aliases: bs_pm,mina,ANY-method bs_pm,mina,character-method

### ** Examples

maize <- new("mina", tab = maize_asv2, des = maize_des2)
maize <- norm_tab(maize, method = "raref")
maize <- fit_tabs(maize)
maize <- bs_pm(maize, group = "Compartment", per = 0.5)



cleanEx()
nameEx("bs_pm")
### * bs_pm

flush(stderr()); flush(stdout())

### Name: bs_pm
### Title: Inferring the network of different group of samples and test
###   significance by permutation.
### Aliases: bs_pm

### ** Examples

maize <- new("mina", tab = maize_asv2, des = maize_des2)
maize <- norm_tab(maize, method = "raref")
maize <- fit_tabs(maize)
maize <- bs_pm(maize, group = "Compartment", per = 0.5)



cleanEx()
nameEx("check_mina")
### * check_mina

flush(stderr()); flush(stdout())

### Name: check_mina
### Title: Check the tab and des file. Return TRUE if they are NULL or both
###   quantitative and descriptive files of same samples are included (i.e.
###   the object is valid).
### Aliases: check_mina
### Keywords: internal

### ** Examples

## Not run: 
##D data(maize)
##D check_mina(maize)
## End(Not run)



cleanEx()
nameEx("check_mina_de")
### * check_mina_de

flush(stderr()); flush(stdout())

### Name: check_mina_de
### Title: Check the object and return TRUE if the object includes
###   descriptive table contains the same samples as quantitative table.
### Aliases: check_mina_de

### ** Examples

## Not run: 
##D data(maize)
##D check_mina_de(maize)
## End(Not run)



cleanEx()
nameEx("check_mina_qu")
### * check_mina_qu

flush(stderr()); flush(stdout())

### Name: check_mina_qu
### Title: Check the object and return TRUE if the object includes
###   quantitative table.
### Aliases: check_mina_qu
### Keywords: internal

### ** Examples

## Not run: 
##D data(maize)
##D check_mina_qu(maize)
## End(Not run)



cleanEx()
nameEx("com_dis-matrix")
### * com_dis-matrix

flush(stderr()); flush(stdout())

### Name: com_dis,matrix,ANY-method
### Title: Calculate the community dissimilarity / distance matrix of the
###   input matrix.
### Aliases: com_dis,matrix,ANY-method com_dis,matrix,character-method

### ** Examples

asv_norm <- norm_tab(maize_asv2, method = "raref", depth = 1000)
asv_dis <- com_dis(asv_norm, method = "bray")



cleanEx()
nameEx("com_dis-mina")
### * com_dis-mina

flush(stderr()); flush(stdout())

### Name: com_dis,mina,ANY-method
### Title: Calculate the community dissimilarity / distance matrix of @norm
###   with 'mina' class object as input.
### Aliases: com_dis,mina,ANY-method com_dis,mina,character-method

### ** Examples

maize <- new("mina", tab = maize_asv2, des = maize_des2)
maize <- norm_tab(maize, method = "total")
maize <- fit_tabs(maize)
maize <- com_dis(maize, method = "bray")



cleanEx()
nameEx("com_dis")
### * com_dis

flush(stderr()); flush(stdout())

### Name: com_dis
### Title: Calculate the community dissimilarity / distance matrix of
###   @norm.
### Aliases: com_dis

### ** Examples

asv_norm <- norm_tab(maize_asv2, method = "raref", depth = 1000)
asv_dis <- com_dis(asv_norm, method = "bray")



cleanEx()
nameEx("com_dis_list")
### * com_dis_list

flush(stderr()); flush(stdout())

### Name: com_dis_list
### Title: List of dissimilarity / distance supported in 'com_dis'.
###   Dissimilarity / distance should be specified by exact string match.
### Aliases: com_dis_list
### Keywords: datasets

### ** Examples

? com_dis_list



cleanEx()
nameEx("com_plot-mina")
### * com_plot-mina

flush(stderr()); flush(stdout())

### Name: com_plot,mina,ANY,ANY,ANY,ANY-method
### Title: Visulization of components distance / dissimilarity in k
###   dimension.
### Aliases: com_plot,mina,ANY,ANY,ANY,ANY-method
###   com_plot,mina,character,ANY,ANY,ANY-method
###   com_plot,mina,character,ANY,ANY,character-method

### ** Examples

maize <- new("mina", tab = maize_asv, des = maize_des)
maize <- norm_tab(maize, method = "raref")
maize <- fit_tabs(maize)
maize <- com_dis(maize, method = "bray")
maize <- dmr(maize)
p1a <- com_plot(maize, match = "Sample_ID", color = "Compartment")
p1b <- com_plot(maize, match = "Sample_ID", d1 = 3, d2 = 4,
color = "Compartment")
p2a <- com_plot(maize, match = "Sample_ID", color = "Host_genotype")
p2b <- com_plot(maize, match = "Sample_ID", d1 = 1, d2 = 3, color =
"Host_genotype")
p3a <- com_plot(maize, match = "Sample_ID", color = "Compartment", shape =
"Soil")
p3b <- com_plot(maize, match = "Sample_ID", d1 = 1, d2 = 4, color =
"Compartment", shape = "Soil")



cleanEx()
nameEx("com_plot")
### * com_plot

flush(stderr()); flush(stdout())

### Name: com_plot
### Title: Visulization of components distance / dissimilarity in k
###   dimension.
### Aliases: com_plot

### ** Examples

maize <- new("mina", tab = maize_asv2, des = maize_des2)
maize <- norm_tab(maize, method = "raref", depth = 5000)
maize <- fit_tabs(maize)
maize <- com_dis(maize, method = "bray")
maize <- dmr(maize)
p1a <- com_plot(maize, match = "Sample_ID", color = "Compartment")
p1b <- com_plot(maize, match = "Sample_ID", d1 = 3, d2 = 4,
color = "Compartment")
p2a <- com_plot(maize, match = "Sample_ID", color = "Host_genotype")
p2b <- com_plot(maize, match = "Sample_ID", d1 = 1, d2 = 3, color =
"Host_genotype")
p3a <- com_plot(maize, match = "Sample_ID", color = "Compartment", shape =
"Soil")
p3b <- com_plot(maize, match = "Sample_ID", d1 = 1, d2 = 4, color =
"Compartment", shape = "Soil")



cleanEx()
nameEx("com_r2-mina")
### * com_r2-mina

flush(stderr()); flush(stdout())

### Name: com_r2,mina,ANY-method
### Title: Function for unexplained variance ratio calculation indicated in
###   Anderson, M.J. 2001. A new method for non-parametric multivariate
###   analysis of variance. Austral Ecology, 26: 32-46.
### Aliases: com_r2,mina,ANY-method com_r2,mina,character-method

### ** Examples

maize <- new("mina", tab = maize_asv2, des = maize_des2)
maize <- norm_tab(maize, method = "raref")
maize <- fit_tabs(maize)
maize <- com_dis(maize, method = "bray")
com_r2(maize, group = c("Compartment", "Soil", "Host_genotype"))



cleanEx()
nameEx("com_r2")
### * com_r2

flush(stderr()); flush(stdout())

### Name: com_r2
### Title: Calculate the unexplained variance ratio using formula indicated
###   in: Anderson, M.J. 2001. A new method for non-parametric multivariate
###   analysis of variance. Austral Ecology, 26: 32-46.
### Aliases: com_r2

### ** Examples

data(maize)
maize <- norm_tab(maize, method = "raref", depth = 5000)
maize <- fit_tabs(maize)
maize <- com_dis(maize, method = "bray")
com_r2(maize, group = c("Compartment", "Soil", "Host_genotype"))



cleanEx()
nameEx("data-hmp")
### * data-hmp

flush(stderr()); flush(stdout())

### Name: data-hmp
### Title: Internal testing data of HMP project, including quantitative
###   table (hmp_otu) and descriptive table (hmp_des) for testing.
### Aliases: data-hmp hmp
### Keywords: data

### ** Examples

data(hmp)



cleanEx()
nameEx("data-maize")
### * data-maize

flush(stderr()); flush(stdout())

### Name: data-maize
### Title: Internal testing data of maize project, vegetative stage samples
###   only, including quantitative table (maize_asv.rds) and descriptive
###   table (maize_des.txt) for testing.
### Aliases: data-maize maize
### Keywords: data

### ** Examples

data(maize)



cleanEx()
nameEx("dmr-matrix")
### * dmr-matrix

flush(stderr()); flush(stdout())

### Name: dmr,matrix-method
### Title: Dimensionality reduction of the distance matrix.
### Aliases: dmr,matrix-method

### ** Examples

## Not run: 
##D data(maize)
##D maize <- norm_tab(maize, method = "raref")
##D maize <- fit_tabs(maize)
##D maize <- com_dis(maize, method = "bray")
##D asv_dis <- maize@dis
##D asv_dis_dmr <- dmr(asv_dis)
##D asv_dis_dmr <- dmr(asv_dis, k = 4)
## End(Not run)



cleanEx()
nameEx("dmr-mina")
### * dmr-mina

flush(stderr()); flush(stdout())

### Name: dmr,mina-method
### Title: Dimensionality reduction of the @dis included in mina.
### Aliases: dmr,mina-method

### ** Examples

## Not run: 
##D data(maize)
##D maize <- norm_tab(maize, method = "raref")
##D maize <- fit_tabs(maize)
##D maize <- com_dis(maize, method = "bray")
##D maize <- dmr(maize)
##D maize <- dmr(maize, k = 4)
## End(Not run)



cleanEx()
nameEx("dmr")
### * dmr

flush(stderr()); flush(stdout())

### Name: dmr
### Title: Dimensionality reduction of community dissimilarity / distance
###   for visulization.
### Aliases: dmr

### ** Examples

maize <- new("mina", tab = maize_asv2, des = maize_des2)
maize <- norm_tab(maize, method = "raref")
maize <- fit_tabs(maize)
maize <- com_dis(maize, method = "bray")
maize <- dmr(maize)



cleanEx()
nameEx("fit_tabs-mina")
### * fit_tabs-mina

flush(stderr()); flush(stdout())

### Name: fit_tabs,mina-method
### Title: Filter the quantitative and descriptive table to make them have
###   the same samples, samples present in both tables are remained. If
###   @norm table exist in the 'mina' object, descriptive table will be
###   filtered again to only keep samples present in @norm.
### Aliases: fit_tabs,mina-method

### ** Examples

{
data(maize)
maize <- fit_tabs(maize)
maize <- norm_tab(maize, method = "raref")
maize <- fit_tabs(maize)
}



cleanEx()
nameEx("fit_tabs")
### * fit_tabs

flush(stderr()); flush(stdout())

### Name: fit_tabs
### Title: Filter the quantitative and descriptive table to make them have
###   the same samples, the intersect samples will be remained.
### Aliases: fit_tabs

### ** Examples

data(maize)
maize <- fit_tabs(maize)
maize <- norm_tab(maize, method = "raref")
maize <- fit_tabs(maize)



cleanEx()
nameEx("get_net_cls_tab-matrix-data.frame-method")
### * get_net_cls_tab-matrix-data.frame-method

flush(stderr()); flush(stdout())

### Name: get_net_cls_tab,matrix,data.frame-method
### Title: Get the cluster table @cls_tab from quantitative table @norm and
###   network clustering results @cls.
### Aliases: get_net_cls_tab,matrix,data.frame-method

### ** Examples

maize <- new("mina", tab = maize_asv2, des = maize_des2)
maize <- norm_tab(maize, method = "raref", depth = 1000)
maize <- fit_tabs(maize)
maize_norm <- maize@norm
maize_adj <- adj(maize_norm, method = "spearman")
maize_cls <- net_cls(maize_adj, method = "ap", cutoff = 0.5)
maize_cls_tab <- get_net_cls_tab(maize_norm, maize_cls)



cleanEx()
nameEx("get_net_cls_tab")
### * get_net_cls_tab

flush(stderr()); flush(stdout())

### Name: get_net_cls_tab
### Title: Get the cluster table @cls_tab from quantitative table @norm and
###   network clustering results @cls.
### Aliases: get_net_cls_tab

### ** Examples

maize <- new("mina", tab = maize_asv2, des = maize_des2)
maize <- norm_tab(maize, method = "raref", depth = 1000)
maize <- fit_tabs(maize)
maize_norm <- maize@norm
maize_adj <- adj(maize_norm, method = "spearman")
maize_cls <- net_cls(maize_adj, method = "ap", cutoff = 0.5)
maize_cls_tab <- get_net_cls_tab(maize_norm, maize_cls)



cleanEx()
nameEx("get_r2-mat")
### * get_r2-mat

flush(stderr()); flush(stdout())

### Name: get_r2,matrix,ANY,ANY-method
### Title: Function for unexplained variance ratio calculation indicated in
###   Anderson, M.J. 2001. A new method for non-parametric multivariate
###   analysis of variance. Austral Ecology, 26: 32-46.
### Aliases: get_r2,matrix,ANY,ANY-method
###   get_r2,matrix,data.frame,ANY-method
###   get_r2,matrix,data.frame,character-method

### ** Examples

maize <- new("mina", tab = maize_asv2, des = maize_des2)
maize <- norm_tab(maize, method = "raref")
maize <- fit_tabs(maize)
maize <- com_dis(maize, method = "bray")
x <- maize@dis
des <- maize@des
get_r2(x, des, group = c("Compartment", "Soil"))



cleanEx()
nameEx("get_r2")
### * get_r2

flush(stderr()); flush(stdout())

### Name: get_r2
### Title: Same function as 'com_r2' with matrix and corresponding
###   descriptive table as input.
### Aliases: get_r2

### ** Examples

maize <- new("mina", tab = maize_asv2, des = maize_des2)
maize <- norm_tab(maize, method = "raref")
maize <- fit_tabs(maize)
maize <- com_dis(maize, method = "bray")
x <- maize@dis
des <- maize@des
get_r2(x, des, group = c("Compartment", "Soil"))



cleanEx()
nameEx("hmp_des")
### * hmp_des

flush(stderr()); flush(stdout())

### Name: hmp_des
### Title: Design file for HMP project, including 2711 samples in total.
### Aliases: hmp_des

### ** Examples

data(hmp_des)




cleanEx()
nameEx("hmp_otu")
### * hmp_otu

flush(stderr()); flush(stdout())

### Name: hmp_otu
### Title: OTU table of HMP project, data downloaded from
###   https://www.hmpdacc.org/hmp/HMQCP/
### Aliases: hmp_otu

### ** Examples

data(hmp_otu)




cleanEx()
nameEx("maize_asv")
### * maize_asv

flush(stderr()); flush(stdout())

### Name: maize_asv
### Title: ASV table of maize project, vegetative stage samples only.
### Aliases: maize_asv

### ** Examples

data(maize_asv)




cleanEx()
nameEx("maize_asv2")
### * maize_asv2

flush(stderr()); flush(stdout())

### Name: maize_asv2
### Title: Subset of ASV table of maize project, ASVs appear in less than
###   100 samples were filtered for later analysis.
### Aliases: maize_asv2

### ** Examples

data(maize_asv2)




cleanEx()
nameEx("maize_des")
### * maize_des

flush(stderr()); flush(stdout())

### Name: maize_des
### Title: Design file of maize project, vegetative stage samples only,
###   including 528 samples in total.
### Aliases: maize_des

### ** Examples

data(maize_des)




cleanEx()
nameEx("maize_des2")
### * maize_des2

flush(stderr()); flush(stdout())

### Name: maize_des2
### Title: Subset of design file of maize project, 313 samples are
###   included.
### Aliases: maize_des2

### ** Examples

data(maize_des2)




cleanEx()
nameEx("mina-class")
### * mina-class

flush(stderr()); flush(stdout())

### Name: mina-class
### Title: Class "mina" includes the quantitative table and descriptive
###   table.
### Aliases: mina-class

### ** Examples

maize <- new("mina", tab = maize_asv, des = maize_des)



cleanEx()
nameEx("net_cls-matrix")
### * net_cls-matrix

flush(stderr()); flush(stdout())

### Name: net_cls,matrix,ANY-method
### Title: Network clustering based on the sparsed adjacacency matrix @adj.
### Aliases: net_cls,matrix,ANY-method net_cls,matrix,character-method

### ** Examples

asv_norm <- norm_tab(maize_asv2, method = "raref", depth = 1000)
asv_adj <- adj(asv_norm, method = "spearman")
asv_cls <- net_cls(asv_adj, method = "mcl")



cleanEx()
nameEx("net_cls-mina")
### * net_cls-mina

flush(stderr()); flush(stdout())

### Name: net_cls,mina,ANY-method
### Title: Network clustering based on the sparsed adjacacency matrix @adj.
### Aliases: net_cls,mina,ANY-method net_cls,mina,character-method

### ** Examples

maize <- new("mina", tab = maize_asv2, des = maize_des2)
maize <- norm_tab(maize, method = "raref", depth = 1000)
maize <- fit_tabs(maize)
maize <- adj(maize, method = "spearman")
maize <- net_cls(maize, method = "mcl", cutoff = 0.4, neg = FALSE)
maize <- net_cls(maize, method = "ap", cutoff = 0.4, neg = FALSE)



cleanEx()
nameEx("net_cls")
### * net_cls

flush(stderr()); flush(stdout())

### Name: net_cls
### Title: Network clustering of sparsed adjacacency matrix @adj.
### Aliases: net_cls

### ** Examples

maize <- new("mina", tab = maize_asv2, des = maize_des2)
maize <- norm_tab(maize, method = "raref", depth = 1000)
maize <- fit_tabs(maize)
maize <- adj(maize, method = "spearman")
maize <- net_cls(maize, method = "mcl", cutoff = 0.4, neg = FALSE)



cleanEx()
nameEx("net_cls_tab-mina-method")
### * net_cls_tab-mina-method

flush(stderr()); flush(stdout())

### Name: net_cls_tab,mina-method
### Title: Get the cluster table @cls_tab from quantitative table @norm and
###   network clustering results @cls.
### Aliases: net_cls_tab,mina-method

### ** Examples

maize <- new("mina", tab = maize_asv2, des = maize_des2)
maize <- norm_tab(maize, method = "raref", depth = 1000)
maize <- fit_tabs(maize)
maize <- adj(maize, method = "spearman")
maize <- net_cls(maize, method = "mcl", cutoff = 0.5)
maize <- net_cls_tab(maize)



cleanEx()
nameEx("net_cls_tab")
### * net_cls_tab

flush(stderr()); flush(stdout())

### Name: net_cls_tab
### Title: Get the cluster table @cls_tab from @norm and @cls.
### Aliases: net_cls_tab

### ** Examples

maize <- new("mina", tab = maize_asv2, des = maize_des2)
maize <- norm_tab(maize, method = "raref", depth = 1000)
maize <- fit_tabs(maize)
maize <- adj(maize, method = "spearman")
maize <- net_cls(maize, method = "ap", cutoff = 0.5)
maize <- net_cls_tab(maize)



cleanEx()
nameEx("net_dis-mina")
### * net_dis-mina

flush(stderr()); flush(stdout())

### Name: net_dis,mina,ANY-method
### Title: Calculate the network distance of @multi and test the
###   significance when @perm is defined.
### Aliases: net_dis,mina,ANY-method net_dis,mina,character-method

### ** Examples

maize <- new("mina", tab = maize_asv2, des = maize_des2)
maize <- norm_tab(maize, method = "raref")
maize <- fit_tabs(maize)
maize <- bs_pm(maize, group = "Compartment")
maize <- net_dis(maize, method = "Jaccard")



cleanEx()
nameEx("net_dis")
### * net_dis

flush(stderr()); flush(stdout())

### Name: net_dis
### Title: Calculate the network distance of @multi and test the
###   significance when @perm is defined.
### Aliases: net_dis

### ** Examples

maize <- new("mina", tab = maize_asv2, des = maize_des2)
maize <- norm_tab(maize, method = "raref")
maize <- fit_tabs(maize)
maize <- bs_pm(maize, group = "Compartment")
## Not run: 
##D maize <- net_dis(maize, method = "spectra", evk = 30)
## End(Not run)



cleanEx()
nameEx("net_dis_indi")
### * net_dis_indi

flush(stderr()); flush(stdout())

### Name: net_dis_indi
### Title: Calculate the network distance of bootstrap and permutation when
###   appliable.
### Aliases: net_dis_indi net_dis_indi,character,ANY-method
###   net_dis_indi,character,character-method

### ** Examples

## Not run: 
##D data(maize)
##D maize <- norm_tab(maize, method = "raref")
##D maize <- fit_tabs(maize)
##D maize <- bs_pm(maize, group = "Compartment", individual = TRUE, out_dir =
##D "./individual_bs_pm/")
##D maize_stat1 <- net_dis_indi(x = "./individual_bs_pm/", method = "spectra")
##D maize_stat2 <- net_dis_indi(x = "./individual_bs_pm/", method = "Jaccard")
##D maize_stat3 <- net_dis_indi(x = "./individual_bs_pm/", method = "spectra",
##D evk = 100, skip = TRUE)
## End(Not run)
## Not run: 
##D data(maize)
##D maize <- norm_tab(maize, method = "raref")
##D maize <- fit_tabs(maize)
##D maize <- bs_pm(maize, group = "Compartment", individual = TRUE, out_dir =
##D "./individual_bs_pm/")
##D maize_stat1 <- net_dis_indi("./individual_bs_pm/", method = "spectra")
##D maize_stat2 <- net_dis_indi("./individual_bs_pm/", method = "Jaccard")
##D maize_stat3 <- net_dis_indi("./individual_bs_pm/", method = "spectra",
##D evk = 100, skip = TRUE)
## End(Not run)



cleanEx()
nameEx("norm_tab-matrix")
### * norm_tab-matrix

flush(stderr()); flush(stdout())

### Name: norm_tab,matrix,character-method
### Title: Normalize the quantitative matrix.
### Aliases: norm_tab,matrix,character-method

### ** Examples

data(maize_asv2)
maize_asv_norm <- norm_tab(maize_asv2, method = "total")
maize_asv_norm <- norm_tab(maize_asv2, method = "raref", depth = 1000,
replace = TRUE, multi = 3)



cleanEx()
nameEx("norm_tab-mina")
### * norm_tab-mina

flush(stderr()); flush(stdout())

### Name: norm_tab,mina,ANY-method
### Title: Normalize the quantitative table with mina input.
### Aliases: norm_tab,mina,ANY-method norm_tab,mina,character-method

### ** Examples

maize <- new("mina", tab = maize_asv2, des = maize_des2)
maize <- norm_tab(maize, method = "raref", depth = 1000, replace = TRUE,
multi = 3)



cleanEx()
nameEx("norm_tab")
### * norm_tab

flush(stderr()); flush(stdout())

### Name: norm_tab
### Title: Normalize the @tab and obtain @norm for later analysis.
### Aliases: norm_tab

### ** Examples

maize <- new("mina", tab = maize_asv2, des = maize_des2)
maize <- norm_tab(maize, method = "total")
maize <- norm_tab(maize, method = "raref")
maize <- norm_tab(maize, method = "raref", depth = 1000, replace = TRUE)
maize <- norm_tab(maize, method = "raref", depth = 1000, replace = TRUE,
multi = 99)



cleanEx()
nameEx("norm_tab_method_list")
### * norm_tab_method_list

flush(stderr()); flush(stdout())

### Name: norm_tab_method_list
### Title: List of normalization methods supported in 'norm_tab'
### Aliases: norm_tab_method_list
### Keywords: datasets

### ** Examples

? norm_tab_method_list



cleanEx()
nameEx("pcoa_plot-list-data.frame-character-ANY-ANY-character-method")
### * pcoa_plot-list-data.frame-character-ANY-ANY-character-method

flush(stderr()); flush(stdout())

### Name: pcoa_plot,list,data.frame,character,ANY,ANY,character-method
### Title: Visulization of components distance / dissimilarity in k
###   dimension.
### Aliases: pcoa_plot,list,data.frame,character,ANY,ANY,character-method

### ** Examples

maize <- new("mina", tab = maize_asv2, des = maize_des2)
maize <- norm_tab(maize, method = "raref")
maize <- fit_tabs(maize)
maize <- com_dis(maize, method = "bray")
maize <- dmr(maize)
asv_dmr <- maize@dmr
des <- maize@des
p1a <- pcoa_plot(asv_dmr, des, match = "Sample_ID", color = "Compartment")
p1b <- pcoa_plot(asv_dmr, des, match = "Sample_ID", d1 = 3, d2 = 4, color =
"Compartment")
p2a <- pcoa_plot(asv_dmr, des, match = "Sample_ID", color = "Host_genotype")
p2b <- pcoa_plot(asv_dmr, des, match = "Sample_ID", d1 = 1, d2 = 3, color =
"Host_genotype")
p3a <- pcoa_plot(asv_dmr, des, match = "Sample_ID", color = "Compartment",
shape = "Soil")
p3b <- pcoa_plot(asv_dmr, des, match = "Sample_ID", d1 = 1, d2 = 4, color =
"Compartment", shape = "Soil")



cleanEx()
nameEx("pcoa_plot")
### * pcoa_plot

flush(stderr()); flush(stdout())

### Name: pcoa_plot
### Title: Visulization of components distance / dissimilarity in k
###   dimension.
### Aliases: pcoa_plot

### ** Examples

maize <- new("mina", tab = maize_asv2, des = maize_des2)
maize <- norm_tab(maize, method = "raref")
maize <- fit_tabs(maize)
maize <- com_dis(maize, method = "bray")
maize <- dmr(maize)
asv_dmr <- maize@dmr
des <- maize@des
p1a <- pcoa_plot(asv_dmr, des, match = "Sample_ID", color = "Compartment")
p1b <- pcoa_plot(asv_dmr, des, match = "Sample_ID", d1 = 3, d2 = 4, color =
"Compartment")
p2a <- pcoa_plot(asv_dmr, des, match = "Sample_ID", color = "Host_genotype")
p2b <- pcoa_plot(asv_dmr, des, match = "Sample_ID", d1 = 1, d2 = 3, color =
"Host_genotype")
p3a <- pcoa_plot(asv_dmr, des, match = "Sample_ID", color = "Compartment",
shape = "Soil")
p3b <- pcoa_plot(asv_dmr, des, match = "Sample_ID", d1 = 1, d2 = 4, color =
"Compartment", shape = "Soil")



cleanEx()
nameEx("sim_par")
### * sim_par

flush(stderr()); flush(stdout())

### Name: sim_par
### Title: Function for community similarity calculation used by 'tina',
###   modified from
###   https://github.com/defleury/Schmidt_et_al_2016_community_similarity/blob/
###   master/functions.community_similarity.R
### Aliases: sim_par
### Keywords: internal

### ** Examples

## Not run: 
##D data(maize_asv)
##D maize_tab <- maize_asv[1 : 1000, 1 : 200]
##D asv <- norm_tab(maize_tab, method = "raref", depth = 100)
##D asv[is.na(asv)] <- 0
##D asv_sparcc <- sparcc(asv, threads = 8, nblocks = 40)
##D tmp.S <- adj(asv_sparcc, method = "spearman")
##D y <- 0.5 * (tmp.S + 1)
##D s <- sim_par(asv_sparcc, y, sim_method = "w_ja", threads = 8, nblocks = 40)
## End(Not run)



cleanEx()
nameEx("sparcc")
### * sparcc

flush(stderr()); flush(stdout())

### Name: sparcc
### Title: Function for 'sparcc' correlation calculation. Modified from
###   Schmidt et al. 2016, find the scripts here and the SparCC paper here.
### Aliases: sparcc
### Keywords: internal

### ** Examples

## Not run: 
##D asv_sparcc <- sparcc(maize_asv2, threads = 2, nblocks = 40)
## End(Not run)



cleanEx()
nameEx("tina-matrix-character-character-method")
### * tina-matrix-character-character-method

flush(stderr()); flush(stdout())

### Name: tina,matrix,character,character-method
### Title: Function for 'tina' dissimilarity calculation. Modified from
###   Schmidt et al., 2016. Person and Spearman could be used for
###   correlation and weighted and unweighted Jaccard could be used for
###   similarity calculation.
### Aliases: tina,matrix,character,character-method

### ** Examples

## Not run: 
##D asv_norm <- norm_tab(maize_asv2, method = "raref", depth = 1000)
##D asv_dis <- com_dis(asv_norm, method = "bray")
##D asv_dis <- com_dis(asv_norm, method = "tina", threads = 8, nblocks = 40)
##D asv_tina <- tina(asv_norm, cor_method = "spearman", sim_method = "w_ja",
##D threads = 8, nblocks = 40)
## End(Not run)



cleanEx()
nameEx("tina")
### * tina

flush(stderr()); flush(stdout())

### Name: tina
### Title: TINA community dissimilarity used in 'com_dis'. Function for
###   'tina' dissimilarity/distance calculation. Modified from Schmidt et
###   al., 2016.
### Aliases: tina

### ** Examples

## Not run: 
##D asv_norm <- norm_tab(maize_asv2, method = "raref", depth = 1000)
##D asv_dis <- com_dis(asv_norm, method = "bray")
##D asv_dis <- com_dis(asv_norm, method = "tina", threads = 8, nblocks = 40)
##D asv_tina <- tina(asv_norm, cor_method = "spearman", sim_method = "w_ja",
##D threads = 8, nblocks = 40)
## End(Not run)



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
