#Set Working Directory
setwd("/proj/vbautchlab/users/arya/zebrafish_Gurung_dataset")

#load needed Libraries
library(Seurat)
library(dplyr)
library(patchwork)
library(SeuratDisk)


#mcherry positive cells should all be endothelial in theory, but clustering and marker analysis is still needed

#read in data from unzipped GEO download 
mcherrypos <- Read10X_h5("data/GSM6138416_filtered_feature_bc_matrix_GFP+mCherry+.h5")

#create the seurat object with basic prefiltering
seurat_object <- CreateSeuratObject(mcherrypos, project = "Gurung_data_set",
                                    min.cells = 3, min.features = 200)

#standard seurat work flow from clustering vignette, 
#some parameters are from the paper, others I determined

#determine the %of mitochondrial DNA
seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")

#look at basic quality control stats before setting cutoffs as needed
VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#cells should have reads from between 91-5456 different genes, with less than 10% of those genes being mitochondrial 
seurat_object <- subset(seurat_object, subset = nFeature_RNA > 96 & nFeature_RNA < 5456 
                        & percent.mt < 10)
#Exclude doublets by setting total counts for a single cell to be less than 61000
seurat_object <- subset(seurat_object, subset = nCount_RNA < 61000)

#see seurat vignette
seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seurat_object)
seurat_object <- ScaleData(seurat_object, features = all.genes)
seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))
#use elbow plot to determine dimesions
ElbowPlot(seurat_object)
seurat_object <- FindNeighbors(seurat_object, dims = 1:12)
seurat_object <- FindClusters(seurat_object, resolution = 0.5)
seurat_object <- RunUMAP(seurat_object, dims = 1:12)
seurat_object <- RunTSNE(seurat_object, dims = 1:12)
DimPlot(seurat_object, reduction = "umap", size = 2, label = TRUE)

#visualize Endothelial Marker expression in the clusters
FeaturePlot(seurat_object, features = cbind("cdh5", "pecam1", "kdr"),
            reduction = "umap")

#pull out clusters with high expression for the given endothelial markers
seurat_object <- subset(seurat_object,
                        idents = cbind("0","1","2","4","5","6","7","11"))


#seurat recommends recalculating Variable Features and rescaling after subsetting if reclustering
seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seurat_object)
seurat_object <- ScaleData(seurat_object, features = all.genes)
seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))
ElbowPlot(seurat_object)
seurat_object <- FindNeighbors(seurat_object, dims = 1:10)
seurat_object <- FindClusters(seurat_object, resolution = 0.5)
seurat_object <- RunUMAP(seurat_object, dims = 1:10)
seurat_object <- RunTSNE(seurat_object, dims = 1:10)
DimPlot(seurat_object, reduction = "umap", size = 2, label = TRUE)

saveRDS(seurat_object, "RDSs/endothelials.rds")
seurat_object <- readRDS("RDSs/endothelials.rds")

markers <- FindAllMarkers(seurat_object)
#only look at markers with adjusted p values less than 0.05
markers <- markers[markers$p_val_adj < 0.05,]
#only look at markers with log 2 fold changes greater tha 0.5 or less than -0.5
markers <- markers[markers$avg_log2FC > 0.5 | markers$avg_log2FC < -0.5,]

#Save Markers for cluster identification
write.csv(markers, file = "markers/endothelial_markers")
