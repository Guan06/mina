library("devtools")
library("dplyr")
load_all("/biodata/dep_psl/grp_rgo/guan/Rpackage_Mina/MINA/")

maize_asv <- readRDS("maize_asv.rds")
maize_des <- read.table("maize_des.txt", header=T, sep="\t")

maize_asv <- as.matrix(maize_asv)
maize_asv <- filter_mat(maize_asv, 100)
maize_asv <- maize_asv[, colSums(maize_asv) > 1000]
maize_des <- maize_des %>% filter(Sample_ID %in% colnames(maize_asv))

maize_asv2 <- maize_asv
maize_des2 <- maize_des

use_data(maize_asv2, overwrite = TRUE)
use_data(maize_des2, overwrite = TRUE)
