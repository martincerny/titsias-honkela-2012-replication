library(tigre)
library(affy)
library(oligo)
library(puma)

source("process.R")

expfiles <- c(paste("array_data/embryo_tc_4_", 1:12, ".CEL", sep=""),
paste("array_data/embryo_tc_6_", 1:12, ".CEL", sep=""),
paste("array_data/embryo_tc_8_", 1:12, ".CEL", sep=""))
expdata <- read.celfiles(expfiles)
pData(expdata) <- data.frame("time.h" = rep(1:12, 3),row.names=rownames(pData(expdata)))


drosophila_mmgmos_exprs <- mmgmos(expdata)
drosophila_gpsim_processed <- processData(drosophila_mmgmos_exprs, experiments=rep(1:3, each=12))
drosophila_gpsim_normalized = normalizeGenes(drosophila_gpsim_processed)

#Reformat the data to separate replicates
yReplicates = array(-1, c(length(drosophila_gpsim_normalized$experiments), length(drosophila_gpsim_normalized$genes), length(drosophila_gpsim_normalized$times)));
yvarReplicates = array(-1, c(length(drosophila_gpsim_normalized$experiments), length(drosophila_gpsim_normalized$genes), length(drosophila_gpsim_normalized$times)));
for(replicate in 1:length(drosophila_gpsim_normalized$experiments)) {
  startIndex = (replicate - 1) * length(drosophila_gpsim_normalized$times) + 1;
  endIndex = replicate * length(drosophila_gpsim_normalized$times);
  yReplicates[replicate,,] = drosophila_gpsim_normalized$y[,startIndex:endIndex];
  yvarReplicates[replicate,,] = drosophila_gpsim_normalized$yvar[,startIndex:endIndex];
}
drosophila_gpsim_normalized$y = yReplicates
drosophila_gpsim_normalized$yvar = yvarReplicates

#Cleanup intermediary data
rm(drosophila_mmgmos_exprs)
rm(expdata)
rm(expfiles)
rm(drosophila_gpsim_processed)
rm(yReplicates)
rm(yvarReplicates)


mapping = loadMapping();

final_training_genes_ids = c(97,216,251,570,577,1325,2732,2733,2735,2773,3430,3888,3900,4133,4394,4512,4654,4795,10433,33509,30900,31313,38134,39039,40089)
#The gene FBgn0026403 (from supplement) is actually coded in the mapping as FBgn0033509 so replaced in the line above

final_training_genes = sprintf("FBgn%07d",final_training_genes_ids)

rm(final_training_genes_ids)

final_training_spots = unlist(lapply(final_training_genes, spotsFromFBgn, mapping=mapping))

twiSpotNames = spotsFromSymbol(mapping, "twi")
allTFSpots =  unlist(lapply(c("tin","bin","twi","bap","Mef2"), function(x) { spotsFromSymbol(mapping,x) } ))

