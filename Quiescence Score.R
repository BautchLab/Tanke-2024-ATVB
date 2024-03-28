#QUIESCENCE SCORE on scRNA-Seq HUVEC DATASET

#load the needed libraries
library(Seurat)
library(dplyr)
library(patchwork)

#Read in Seurat Object of the HUVEC Static Flow Dataset(GSE151867)
stat <- readRDS("HUVEC_v3.RDS")

#change to Epithelial CSVs as needed
Quiupgenes <- read.csv("Endothelial Up.csv", header = 1)$Gene
Quidwgenes <- read.csv("Endothelial Down.csv", header = 1)$Gene

#Use @scale.data to give equal weight to all quiescent genes regardless of gene expression level
calscore <- function (EC){
  df_up <- EC@assays$RNA@scale.data[rownames(EC) %in% Quiupgenes,]
  df_dw <- EC@assays$RNA@scale.data[rownames(EC) %in% Quidwgenes,]
  score_up <- apply(df_up,2,function(x) mean(x, na.rm = T))
  EC@meta.data$Qui_score_up <- score_up
  score_dw <- apply(df_dw,2,function(x) mean(x, na.rm = T))
  score <- as.data.frame(cbind(score_up,score_dw, (score_up - score_dw)))
  colnames(score) <- c("Qui_score_up","Qui_score_dw","Qui_score")
  EC@meta.data$Qui_score_up <- score$Qui_score_up
  EC@meta.data$Qui_score_dw <- score$Qui_score_dw
  EC@meta.data$Qui_score <- score$Qui_score
  return (EC)
}

# now calculate quiescence score for the dataset using the function "calscore" defined above
stat <- calscore(stat)

# a Seurat Object with the quiescence scores as a meta.data column should be returned
#Take a quick look
scores <- c("Qui_score_up","Qui_score_dw","Qui_score")
VlnPlot(stat, features = scores, ncol= 3, pt.size= 0)

#write out Quiescence Score for Plotting
stat[["subpopulation"]]<-Idents(stat)
temp<-split(stat@meta.data$Qui_score,stat@meta.data$subpopulation)
indx <- sapply(temp, length)
scores <- do.call(cbind,lapply(temp, `length<-`, max(indx)))
write.csv(scores,"Qui_score_by_subpopulations.csv")