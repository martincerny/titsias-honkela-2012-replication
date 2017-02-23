library(drosgenome1.db)
library(tigre)
library(affy)
library(oligo)
library(puma)


expfiles <- c(paste("array_data/embryo_tc_4_", 1:12, ".CEL", sep=""),
paste("array_data/embryo_tc_6_", 1:12, ".CEL", sep=""),
paste("array_data/embryo_tc_8_", 1:12, ".CEL", sep=""))
expdata <- read.celfiles(expfiles)
pData(expdata) <- data.frame("time.h" = rep(1:12, 3),row.names=rownames(pData(expdata)))


drosophila_mmgmos_exprs <- mmgmos(expdata)
drosophila_mmgmos_fragment <- drosophila_mmgmos_exprs
drosophila_gpsim_fragment <- processData(drosophila_mmgmos_fragment, experiments=rep(1:3, each=12))


final_training_genes_ids = c(97,216,251,570,577,1325,2732,2733,2735,2773,3430,3888,3900,4133,4394,4512,4654,4795,10433,26403,30900,31313,38134,39039,40089)

final_trainining_genes = sprintf("FBgn%07d",final_training_genes_ids)

#Mapping to recognized synonym:
#570 -> 260400
mapToSynonyms <- function(x)
{
  x[x == "FBgn0000570"] <- "FBgn0260400"
}

final_training_spots = mget(mapToSynonyms(final_trainining_genes), env=drosgenome1FLYBASE2PROBE)