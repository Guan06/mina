library("devtools")
library("dplyr")
load_all("/biodata/dep_psl/grp_rgo/guan/Rpackage_Mina/MINA/")

hmp_otu <- readRDS("hmp_otu.rds")
hmp_otu <- as.matrix(hmp_otu)

hmp_des <- read.table("hmp_des.txt", header=T, sep="\t")
hmp_des <- hmp_des %>% filter(Sample_ID %in% colnames(hmp_otu))

hmp_otu <- filter_mat(hmp_otu, 10)
hmp_otu <- hmp_otu[, colSums(hmp_otu) > 500]
hmp_des <- hmp_des %>% filter(Sample_ID %in% colnames(hmp_otu))

use_data(hmp_otu, overwrite = TRUE)
use_data(hmp_des, overwrite = TRUE)

load_all("/biodata/dep_psl/grp_rgo/guan/Rpackage_Mina/MINA/")
hmp <- new("mina", tab = hmp_otu, des = hmp_des)
use_data(hmp, overwrite = TRUE)
