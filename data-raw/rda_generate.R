library("devtools")
library("dplyr")
maize_asv <- readRDS("maize_asv.rds")
maize_asv <- as.matrix(maize_asv)

maize_des <- read.table("maize_des.txt", header=T, sep="\t")
maize_des <- maize_des %>% filter(Sample_ID %in% colnames(maize_asv))

use_data(maize_asv, overwrite = TRUE)
use_data(maize_des, overwrite = TRUE)
