library("devtools")
library("dplyr")
maize_asv <- readRDS("maize_asv.rds")
maize_asv <- as.matrix(maize_asv)

maize_des <- read.table("maize_des.txt", header=T, sep="\t")
maize_des <- maize_des %>% filter(Sample_ID %in% colnames(maize_asv))

maize_asv <- filter_mat(maize_asv, 10)
maize_asv <- maize_asv[, colSums(maize_asv) > 500]
maize_des <- maize_des %>% filter(Sample_ID %in% colnames(maize_asv))

use_data(maize_asv, overwrite = TRUE)
use_data(maize_des, overwrite = TRUE)

load_all("/biodata/dep_psl/grp_rgo/guan/Rpackage_Mina/MINA/")
maize <- new("mina", tab = maize_asv, des = maize_des)
use_data(maize, overwrite = TRUE)
