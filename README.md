---
title: "Rat snRNA-seq data analysis"
format: html
editor: visual
R version: 4.4.1 (2024-16-14) -- "Race for Your Life"
---

### Rat snRNA-seq data analysis 

Systemic and breast chronic inflammation and hormone disposition promote a tumor-permissive locale for breast cancer in older women

![image](https://github.com/user-attachments/assets/a5a4ad70-20e1-4c36-976d-5cfe9c477cc4)

### Necessary input data files. (See Preparation2 below)

-   snRNA-seq processed count data 
-   Age group annotation data

## Preparation1. Install and load the necessary R packages.
```{r}
library(Rcpp)
library(BiocManager) # CRAN # install.packages("BiocManager")
library(Seurat, warn.conflicts = FALSE)  # install.packages("Seurat")  Seurat version 5.1.0 https://satijalab.org/seurat/
library(SeuratObject, warn.conflicts = FALSE) 
library(stringr) # CRAN
library(ggplot2) # CRAN
library(dplyr, warn.conflicts = FALSE) # CRAN
library(tidyr) # CRAN
library(data.table, warn.conflicts = FALSE) # CRAN
library(dittoSeq)  # Bioconductor - it requires "nloptr" CRAN package # BiocManager::install("dittoSeq")
library(glmGamPoi, warn.conflicts = FALSE) # Bioconuctor package for SCTransform.  
# library(cowplot) # CRAN, for plot_grid()

library(patchwork)  # CRAN. to use plot_annotation()
library(harmony)  # CRAN.   install.packages("harmony") to use RunHarmony() https://github.com/immunogenomics/harmony
```

## Preparation2. Set working directory and direct your input data files
```{r}
# Your working directory that this code file
dir <- dirname(rstudioapi::getSourceEditorContext()$path); 
setwd(dir); print(dir)

## You can download this file from NCBI GEO GSEOOOOOOOOO
FilteredFeatureCountDir <- "/Users/sanghoonlee/Library/CloudStorage/OneDrive-UniversityofPittsburgh/H08_AgeStudy_Rat_snRNAseq_Neil/04_CellRangerCount"

## You can download this file from NCBI GEO GSEOOOOOOOOO
AgegroupAnnotFile <-  "/Users/sanghoonlee/Library/CloudStorage/OneDrive-UniversityofPittsburgh/H08_AgeStudy_Rat_snRNAseq_Neil/04_CellRangerCount/Rat_scRNAseq_AgeGroup.txt"
```

## Section 1. Read 10X Filtered feature Count matrix andset up a Seurat object. 
```{r}
AgingRat_ERp_Dir <- list.dirs(path=FilteredFeatureCountDir); 
AgingRat_ERp_Dir <- AgingRat_ERp_Dir[-1]; print(AgingRat_ERp_Dir) # remove top dir name. 

SeuratObjList <- list(); ForLoopNumb<-0;
for (EachDir in AgingRat_ERp_Dir) {
      # EachDir <- AgingRat_ERp_Dir[1]
      print(paste0("Each dir: ", EachDir))
      ForLoopNumb <- ForLoopNumb+1;
      
      expression_matrix <- Read10X(data.dir = EachDir, gene.column=2)  # 23098g 17490  # gene read count matrix.
      # Should show exactly the file names of "barcodes.tsv.gz", "features.tsv.gz"(genes.tsv), and "matrix.mtx.gz"   
      # Feature first column will determine the rownames in "seurat_object@assays$RNA@data" whether ENSG ID or GeneSymbol
      # View(expression_matrix) # This has only expression matrix. 
      seurat_object = CreateSeuratObject(counts = expression_matrix)  # 33538   3406 cells
      #dim(seurat_object@assays$RNA@data)   # 33538 genes 3406 cells
      
      SeuratObjName <- gsub("(.*)/", "", EachDir)
      seurat_object$orig.ident <- SeuratObjName
      ################################################################################
      ## Step1. QC and Filter out mitocondrial genes
      ################################################################################
      seurat_object[["percent.mt"]]  <- PercentageFeatureSet(seurat_object, pattern="^MT-")
      seurat_object[["percent.rbp"]]  <- PercentageFeatureSet(seurat_object, pattern="^RP[SL]")   ## 33538  3406 cells
      #VlnPlot(seurat_object, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rbp"), ncol = 4)
      seurat_object <- subset(seurat_object, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA > 400 & percent.mt < 15)
      dim(seurat_object) # 23098 17348
      ################################################################################
      ## Step2. SCTransform
      ################################################################################
      seurat_object <- SCTransform(seurat_object, method="glmGamPoi", vars.to.regress="percent.mt", verbose=FALSE)  #method="glmGamPoi" requires "glmGamPoi" R package
      ################################################################################
      ## Step3. Remove doublet
      ################################################################################
      # seurat_objectClean <- subset(seurat_object, cells=rownames(seurat_object@meta.data)[which(seurat_object@meta.data$DF.classification == "Singlet")])
      # table(seurat_objectClean$DF.classification) # Singlet: 2960
      seurat_objectClean <- seurat_object
      ################################################################################
      ## Step4. Store Seurat objects into a list
      ################################################################################
      #DefaultAssay(seurat_objectClean) <- "RNA" ## I don't need this.
      if(ForLoopNumb==1) {
        MyFirstSeuratObj <- seurat_objectClean                 # This has the first object
      } else {
        SeuratObjList[[SeuratObjName]] <- seurat_objectClean   # This has 19 objectd
      }
      #if(ForLoopNumb==2) break;
}
```

# Section 2. Merge, SCTransform and Harmony integration

## Step2a.  Merge Seurat objects.   # https://satijalab.org/seurat/articles/merge_vignette.html
```{r}
## Merge takes about 10 minutes
Seurat_Merge <- merge(unlist(MyFirstSeuratObj) , unlist(SeuratObjList), add.cell.ids=c(unique(MyFirstSeuratObj$orig.ident),names(SeuratObjList)))   # this works. 
saveRDS(Seurat_Merge, file="SeuratObj_Rat_snRNAseq_AfQC.rds")

```

## Step2b	SCTransform to regress out mitochondrial read percentage per cell  - This requires a lot of memory and takes. You may need to use supercomputer
```{r}
SeuratObject_Rat_snRNAseq_SCTransform <- Seurat::SCTransform(Seurat_Merge, method="glmGamPoi", vars.to.regress="percent.mt", verbose=FALSE)   # This takes about 10 min. # This step requires a lot of memory. It may not run in your local computer. You will need a supercomputer.
saveRDS(SeuratObject_Rat_snRNAseq_SCTransform, "SeuratObject_MergeSCTransform_Rat_snRNAseq.rds")

### If SCTransform fails in you computer, run this line to load Seurat Object.  
# SeuratObject_Rat_snRNAseq <- readRDS("SeuratObject_MergeSCTransform_Rat_snRNAseq.rds")

##### $$$$$ ###### $$$$$$ ###### Supplementary Figure. UMAP plot of integrated datasets by CaseID - Before Harmony Integration ##### $$$$$ ###### $$$$$ #####
SeuratObject_Rat_snRNAseq_SCTransformPCA <- Seurat::RunPCA(SeuratObject_Rat_snRNAseq_SCTransform, npcs = 30, verbose = F)
SeuratObject_Rat_snRNAseq_SCTransformUMAP <- SeuratObject_Rat_snRNAseq_SCTransformPCA %>% Seurat::RunUMAP(reduction="pca", dims=1:30, verbose=F)
MyDimplot_SeuratMergeSCT_ByCaseID_BfHarmony <- Seurat::DimPlot(SeuratObject_Rat_snRNAseq_SCTransformUMAP,reduction="umap",group.by="orig.ident", label=FALSE) + 
  patchwork::plot_annotation(title="Two datasets UMAP before Harmony by CaseID") +
  theme(axis.text.x=element_text(vjust=0.6, size=25,angle=0), axis.text.y=element_text(vjust=0.6, size=25,angle=0))
ggsave(MyDimplot_SeuratMergeSCT_ByCaseID_BfHarmony, height=8,width=11, dpi=300, filename=paste0("FigS1A_OutUMAP_Rat_snRNAseq_ByCaseID_BfHarmony.pdf"), useDingbats=FALSE)
##### ----------- ============== $$$$$$$$$$$$$  ##### ----------- ============== $$$$$$$$$$$$$  ##### ----------- ============== $$$$$$$$$$$$$  
```

## Step2c	Normalization and Scaling, and make DimPlot
```{r}
SeuratObject_Rat_snRNAseq <- Seurat::RunPCA(SeuratObject_Rat_snRNAseq_SCTransform, npcs = 30, verbose = F)

#SeuratObject_Rat_snRNAseq$DatasetID <- ifelse( grepl("CID",SeuratObject_Rat_snRNAseq$orig.ident), "GSE176078",   ifelse( grepl("BIOKEY",SeuratObject_Rat_snRNAseq$orig.ident), "EGA6608", "NotAvail") )  
#table(SeuratObject_Rat_snRNAseq$DatasetID) # GSE176078: 30959 cells    EGA6608: 41120
```

## Step2d.	Integrate using R package Harmony   Don't forget assay.use="SCT"   - Takes 10 minutes
```{r}
SeuratObject_Rat_snRNAseq <- SeuratObject_Rat_snRNAseq %>% harmony::RunHarmony("orig.ident", plot_convergence=T, assay.use="SCT")  
saveRDS(SeuratObject_Rat_snRNAseq, file="SeuratObject_Rat_snRNAseq_Harmony.rds")   

### If Harmony integration fails in you computer, run this line to load Seurat Object.  
# SeuratObject_Rat_snRNAseq <- readRDS("SeuratObject_Rat_snRNAseqIntegration.rds")

##################### Harmony plot After Harmony Integration  #############################
p1 <- Seurat::DimPlot(object=SeuratObject_Rat_snRNAseq, reduction="harmony", pt.size=.1, group.by="orig.ident") + NoLegend()
p2 <- Seurat::VlnPlot(object=SeuratObject_Rat_snRNAseq, features="harmony_1", group.by ="orig.ident", pt.size = .1) + NoLegend()
# plot_grid(p1,p2)
ggsave(p1, height=8,width=8, dpi=300, filename=paste0("Fig1_OutHarmonyplot_Rat_snRNAseq_MergeSCTHarmony.pdf"), useDingbats=FALSE)
```

# Section 3. UMAP reduction, clustering cells, and define cell types.

## Step3a.	Read  clinical data  => add it to Seurat object
```{r}
## Read Agegroup Annot file 
AgegroupAnnot <- data.table::fread(AgegroupAnnotFile, header=TRUE, stringsAsFactors=FALSE);dim(AgegroupAnnot); print(AgegroupAnnot) # 6 2
#                 CaseID AgeGroup
#                 <char>   <char>
# 1: Lee_021924_Nuclei1     Aged
# 2: Lee_021924_Nuclei2     Aged
# 3: Lee_021924_Nuclei3     Aged
# 4: Lee_021924_Nuclei4    Young
# 5: Lee_021924_Nuclei5    Young
# 6: Lee_021924_Nuclei6    Young

## Change active.ident with orig.ident
# SeuratObject_Rat_snRNAseq <- SetIdent(SeuratObject_Rat_snRNAseq, value = SeuratObject_Rat_snRNAseq$orig.ident)
# table(SeuratObject_Rat_snRNAseq@active.ident)
# Lee_021924_Nuclei1 Lee_021924_Nuclei2 Lee_021924_Nuclei3 Lee_021924_Nuclei4 Lee_021924_Nuclei5 Lee_021924_Nuclei6 
# 17348              19899              13545              16216              13157              13964 

## Get metadata of seurat object and inner_join with Agegroup Annot data. 
RatsnRNAseqMetadata <- SeuratObject_Rat_snRNAseq@meta.data %>% data.frame %>% dplyr::mutate(CaseID = orig.ident); dim(RatsnRNAseqMetadata) # 94129    9
table(RatsnRNAseqMetadata$CaseID %in% AgegroupAnnot$CaseID) # All TRUE 94129
# RatsnRNAseqMetadata_AgeGroup <- dplyr::inner_join(RatsnRNAseqMetadata, AgegroupAnnot); dim(RatsnRNAseqMetadata_AgeGroup) # 47066     9  # this doesn't work. I don't know why
RatsnRNAseqMetadata$AgeGroup <- ifelse(RatsnRNAseqMetadata$CaseID %in% c("Lee_021924_Nuclei1","Lee_021924_Nuclei2","Lee_021924_Nuclei3"), "Aged", "Young")
dim(RatsnRNAseqMetadata); table(RatsnRNAseqMetadata$AgeGroup);  # 94129 9  # Aged 50792  Young 43337

## Add metadata to seurat object. 
SeuratObject_Rat_snRNAseq <- AddMetaData(object=SeuratObject_Rat_snRNAseq, metadata=c(RatsnRNAseqMetadata$AgeGroup), col.name=c("AgeGroup"))
table(SeuratObject_Rat_snRNAseq$AgeGroup) #  # Aged 50792  Young 43337
```

```{r}
################################################################################
## Step4b.	UMAP Reduction and clustering
################################################################################
SeuratObject_Rat_snRNAseq <- SeuratObject_Rat_snRNAseq %>% Seurat::RunUMAP(reduction="harmony", dims=1:30, verbose=F) %>% FindNeighbors(reduction="harmony", k.param=15, dim=1:30)  
SeuratObject_Rat_snRNAseq <- SeuratObject_Rat_snRNAseq %>% Seurat::FindClusters(resolution=1.5) %>% identity()
table(SeuratObject_Rat_snRNAseq@active.ident) # Res1.0, 31 clusters; Res1.5, 40 clusters
saveRDS(SeuratObject_Rat_snRNAseq, file="SeuratUMAPCluster_Rat_snRNAseq_Res1.5PC30KP15.rds") 
# SeuratObject_Rat_snRNAseq <- readRDS(file="SeuratUMAPCluster_Rat_snRNAseq_Res1.5PC30KP15.rds")

SeuratObject_Rat_snRNAseq$CaseID <- SeuratObject_Rat_snRNAseq$orig.ident; table(SeuratObject_Rat_snRNAseq$CaseID )
table(SeuratObject_Rat_snRNAseq@active.ident)
# Lee_021924_Nuclei1 Lee_021924_Nuclei2 Lee_021924_Nuclei3 Lee_021924_Nuclei4 Lee_021924_Nuclei5 Lee_021924_Nuclei6 
# 17348              19899              13545              16216              13157              13964 

## Res 1.0
# 0     1     2     3     4     5     6     7     8     9    10    11    12    13    14    15    16    17    18    19    20    21    22    23    24    25    26    27    28    29    30 
# 11925  9473  7378  5589  5559  5173  4610  4576  4014  3752  3291  3151  3125  3108  2802  2329  1950  1931  1787  1610  1332  1281   970   819   803   800   546   257   101    47    40 
## Res 1.5
# 0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   20   21   22   23   24   25   26   27   28   29   30   31   32   33   34   35   36   37 
# 9957 7158 6955 5721 5054 4669 3994 3792 3775 3758 3715 3180 3095 3059 2785 2475 2326 1997 1972 1959 1745 1605 1282 1148 1061  973  820  787  742  662  568  487  257  195   98   94   60   51 
# 38   39 
# 51   47 

### UMAP by cluster ID
UMAP_ByClusterNumber <- Seurat::DimPlot(SeuratObject_Rat_snRNAseq,reduction="umap", label=TRUE, label.size=8) + patchwork::plot_annotation(title="UMAP_Harmony and Clustering number") +
  theme(axis.text.x=element_text(vjust=0.6, size=25,angle=0), axis.text.y=element_text(vjust=0.6, size=25,angle=0))
ggsave(UMAP_ByClusterNumber, height=8,width=11, dpi=300, filename=paste0("UMAPplot_ByClusterNumber_Rat_snRNAseq_MergeSCTHarmony_Res1.5.pdf"), useDingbats=FALSE)

##### ----------- ============== $$$$$$$$$$$$$  ##### ----------- ============== $$$$$$$$$$$$$  ##### ----------- ============== $$$$$$$$$$$$$  
##### $$$$$ ###### $$$$$$ ###### Supplementary Figure 1A. UMAP plot of integrated datasets by CaseID - Checking Harmony Integration ##### $$$$$ ###### $$$$$ #####
MyDimplot_SeuratMergeSCTHarmony_ByCaseID <- Seurat::DimPlot(SeuratObject_Rat_snRNAseq,reduction="umap",group.by="CaseID", label=FALSE) + patchwork::plot_annotation(title="UMAP_Harmony and CaseID") +
  theme(axis.text.x=element_text(vjust=0.6, size=25,angle=0), axis.text.y=element_text(vjust=0.6, size=25,angle=0))
ggsave(MyDimplot_SeuratMergeSCTHarmony_ByCaseID, height=8,width=11, dpi=300, filename=paste0("FigS1A_OutUMAP_Rat_snRNAseq_ByCaseID.pdf"), useDingbats=FALSE)
##### ----------- ============== $$$$$$$$$$$$$  ##### ----------- ============== $$$$$$$$$$$$$  ##### ----------- ============== $$$$$$$$$$$$$  

### UMAP by AgeGroup
UMAP_ByAgeGroup <- Seurat::DimPlot(SeuratObject_Rat_snRNAseq,reduction="umap",group.by="AgeGroup",  label=FALSE, label.size=8) + patchwork::plot_annotation(title="UMAP_Harmony and AgeGroup") +
  theme(axis.text.x=element_text(vjust=0.6, size=25,angle=0), axis.text.y=element_text(vjust=0.6, size=25,angle=0))
ggsave(UMAP_ByAgeGroup, height=8,width=11, dpi=300, filename=paste0("UMAPplot_ByAgeGroup_Rat_snRNAseq_MergeSCTHarmony_Res1.5.pdf"), useDingbats=FALSE)

```

```{r}
################################################################################
## Step4c. Define cell types by marker gene expression   ####
################################################################################
# MyMainMarker <- c("PECAM1","RAMP2","FLT1","CLDN5",   "EPCAM","KRT19","KRT18","CD24",  "PDGFRB","C1R","DCN","COL1A1",   "ACTA2",  "TPSB2","TPSAB1","CPA3",
#                   "CD68","LYZ","TYROBP",   "CD83","MS4A1","MZB1","CD79A",   "CD2","CD3E","CD3D","CD3G","IL7R") # I should include "CD79A" B cell marker. 
## Check existence of gene in single cell data. 
SeuratObject_Rat_snRNAseq[["SCT"]]@data[1:5,1:3]
table("Cd" %in% rownames(SeuratObject_Rat_snRNAseq[["SCT"]]@data)) # 
grep("Il7r",rownames(SeuratObject_Rat_snRNAseq[["SCT"]]@data), value=TRUE) # 

### Main markers 
MyMainMarker <- c("Fap", "Cdh5","Pecam1","Flt1",      "Epcam","Krt18","Krt8","Cd24",  
                  "Tp63",  # Myoepithelial, Tp63 Double negative cell  "Myl9",
                  "Pdgfrb","C1r","Dcn","Col1a1",  # Fibroblast
                  "Acta2",  # Myofibroblast 
                  "Cd34",  # embryonic fibroblast. 
                  
                  "Cd68","Tyrobp","Cd14", # Myeloid
                  "Csf1r","Ace","Adgre1", "Cd47",   # Monocyte
                  "Fcgr3a","Lst1", # CD16+ Monocyte
                  "Apoe","Fabp5","Egr1","Cx3cr1","Slc2a1", # Macrophage 
                  "Msr1","Mrc1",  # M2 Macrophage
                  "Cd86",    # M1 Macrophage
                  "Cd83",  ## Dendritic, Myeloid
                  "Irf7", # plasmacytoid DC 
                  "Flt3","Gzmb","Itgax","Thbd","Fscn1",    # Dendritic cell.  
                  
                  "Ptprc",  # Immune cells or NK cells
                  "Cd2","Cd3e","Cd3d","Cd3g","Il7r", "Ncam1","Klrd1","Il2rb",  # NKT cells. "Fcgr3a" NKT and Macrophage. 
                         
                  "ENSRNOG00000071219","Stat4", #CD4 T cell marker ENSRNOG00000071219 is CD4 
                  "Gzma","Cd8a" # CD8 T cell. 
                  ) # I should include "CD79A" B cell marker. 
        
MyDotPlot_MainType <- Seurat::DotPlot(SeuratObject_Rat_snRNAseq, features=MyMainMarker)+ theme(axis.text.x=element_text(vjust=0.6, size=15,angle=0), axis.text.y=element_text(vjust=0.6, size=15,angle=0)) +
  scale_colour_gradient2(low = "#1515FA", mid = "#FFFAE2", high = "#ff0000") + coord_flip()    # scale_color_viridis_c()    #ffe272
MyDotPlot_MainType

ggsave(MyDotPlot_MainType, height=16,width=17, dpi=300, filename=paste0("OutDotplot_ByMainMarker_HarmoneyRatsnRNAseq_Res1.5_40Cluster_MoreMarker.pdf"), useDingbats=FALSE)


### Epithelial markers  Li et al. 2020 Cell Report
EpithelialMarker <- c("Tp63","Krt17","Krt5","Acta2","Mylk","Myh11",  # Myoepithelial, Tp63 Double negative cell  "Myl9",
                  "Krt18","Krt8",  # Luminal "Krt19" doesn't exist. 
                  "Prlr","Cited1","Pgr","Esr1","Prom1", # Hormone sensitive
                  "Mfge8","Csn3","Wfdc18","Elf5","Ltf", # secretory alveolar  "Trf" doesn't exist. 
                  "Kit","Adh1a3","Cd14", # luminal progenitor
                  "Wap","Glycam1","Olah" # Late milk
                  ) #
MyDotPlot_Epithelial <- Seurat::DotPlot(SeuratObject_Rat_snRNAseq, features=EpithelialMarker)+ theme(axis.text.x=element_text(vjust=0.6, size=15,angle=0), axis.text.y=element_text(vjust=0.6, size=15,angle=0)) +
          scale_colour_gradient2(low = "#1515FA", mid = "#FFFAE2", high = "#ff0000") + coord_flip()    # scale_color_viridis_c()    #ffe272
MyDotPlot_Epithelial

ggsave(MyDotPlot_Epithelial, height=6,width=14, dpi=300, filename=paste0("OutDotplot_ByEpithelialMarker_HarmoneyRatsnRNAseq_Res1.5_40Cluster_MoreMarker.pdf"), useDingbats=FALSE)


### Stromal markers  Li et al. 2020 Cell Report
StromalMarker <- c("Col1a1","Col1a2","Col3a1","Fn1", # Fibroblast
                      "Pecam1","Cdh5","Eng","Sox17","Sele", # Vascular endothelial
                      "Rgs5","Des","Notch3", # Pericyte
                      "Mmrn1","Prox1","Flt4", # Lymphatic endothelial "Ccl21a",  doesn't exist
                      "Ptprc","Cd74","Lyz2","Napsa","Traf1","Flt3", # Dendritic cell
                      "Csf1r","Adgre1","Ms4a7","Mrc1","Cd209f","Cd163", #Macrophage Ma  "Fcgr3", doesn't exist
                      "Mmp12","Mmp13", "Spic", # Macrophage Mb
                      "Cd3g","Cd3d","Cd3e","Gzma","Ncr1","Itgae", "Cd8a","Cd8b", # NK cells
                      "ENSRNOG00000071219","Cd79a","Cd79b")  # ENSRNOG00000071219 is Cd4. "Cd8b1", doesn't exist
                      
MyDotPlot_Stromal <- Seurat::DotPlot(SeuratObject_Rat_snRNAseq, features=StromalMarker)+ theme(axis.text.x=element_text(vjust=0.6, size=15,angle=0), axis.text.y=element_text(vjust=0.6, size=15,angle=0)) +
  scale_colour_gradient2(low = "#1515FA", mid = "#FFFAE2", high = "#ff0000") + coord_flip()    # scale_color_viridis_c()    #ffe272
MyDotPlot_Stromal

ggsave(MyDotPlot_Stromal, height=10,width=15, dpi=300, filename=paste0("OutDotplot_ByStromalMarker_HarmoneyRatsnRNAseq_Res1.5_40Cluster_MoreMarker.pdf"), useDingbats=FALSE)


## FeaturePlot - To pull gene expression cells to front. When features are burried,  it is useful. 
MyFeaturePlot_SingelGene <- Seurat::FeaturePlot(SeuratObject_Rat_snRNAseq, raster=TRUE,features=c("Epcam","Krt18","Fcgr3a", "Cd14","Cd68","Ms4a1","Pecam1","Il7r","Cd3e","Cd8a"), cols=c("lightgrey","lightgrey","darkblue"), order=TRUE, pt.size=1.5) 
ggsave(MyFeaturePlot_SingelGene, height=8,width=14, dpi=300, filename=paste0("OutFeatureUMAP_ByMainMarker_HarmoneyRatsnRNAseq_Res1.5PC30KP15.pdf"), useDingbats=FALSE)

# MyFeaturePlot_CD14 <- Seurat::FeaturePlot(SeuratObject_Rat_snRNAseq, raster=TRUE,features=c("Cd14"), cols=c("lightgrey","blue"), order=TRUE, pt.size=1.8) 
# ggsave(MyFeaturePlot_CD14, height=3.5,width=4, dpi=300, filename=paste0("OutFeatureUMAP_ByMainMarker_HarmoneyRatsnRNAseq_Res1.0PC30KP15_CD14.pdf"), useDingbats=FALSE)


new.cluster.ids_Main<- c("CancerEpithelial","CancerEpithelial", "CancerEpithelial","CancerEpithelial","CancerEpithelial","CancerEpithelial",     "CancerEpithelial","CancerEpithelial","CancerEpithelial","Myoepithelial","CancerEpithelial",
                         "CancerEpithelial","CancerEpithelial","CancerEpithelial","CancerEpithelial","Myeloid",                                "NKTcell","CancerEpithelial","CancerEpithelial", "CancerEpithelial","Myoepithelial",
                         "Myeloid","NKTcell","CancerEpithelial","CancerEpithelial","Fibroblast",                                          "CancerEpithelial","DendriticCell","Myeloid","CancerEpithelial","Endothelial", 
                         "CancerEpithelial","Fibroblast","Endothelial","CancerEpithelial","CancerEpithelial",                                   "CancerEpithelial","DendriticCell","CancerEpithelial","Myoepithelial") # 40 clusters

names(new.cluster.ids_Main) <- levels(SeuratObject_Rat_snRNAseq)
SeuratObject_Rat_snRNAseq <- Seurat::RenameIdents(SeuratObject_Rat_snRNAseq, new.cluster.ids_Main)
table(SeuratObject_Rat_snRNAseq@active.ident); sum(table(SeuratObject_Rat_snRNAseq@active.ident)) # 94,129
# CancerEpithelial    Myoepithelial       ImmuneCell      StromalCell 
#     77318             5550             9268             1993 
# CancerEpithelial    Myoepithelial          Myeloid          NKTcell       Fibroblast    DendriticCell      Endothelial 
#       77318             5550             4822             3608             1230              838              763 

## ==== ##### ====== ##### Add the new cell type definition to Metdata  ## ==== ##### ====== ##### 
SeuratObject_Rat_snRNAseq <- Seurat::AddMetaData(object=SeuratObject_Rat_snRNAseq, metadata=c(SeuratObject_Rat_snRNAseq@active.ident), col.name=c("CellTypeByMarker_RatsnRNAseq"))
saveRDS(SeuratObject_Rat_snRNAseq, file="SeuratObject_Harmony_Rat_snRNAseq_CellTypeAssign.rds")    

##### ----------- ============== $$$$$$$$$$$$$  ##### ----------- ============== $$$$$$$$$$$$$  ##### ----------- ============== $$$$$$$$$$$$$  
##### $$$$$ ###### $$$$$$ ###### Figure 2B. UMAP plot of integrated datasets by cell type definition ##### $$$$$ ###### $$$$$ #####
MyDimplot_CellType <- Seurat::DimPlot(SeuratObject_Rat_snRNAseq, reduction="umap", group.by="CellTypeByMarker_RatsnRNAseq", pt.size=.1, label=TRUE, label.size=8) + NoLegend(); # label=TRUE,
ggsave(MyDimplot_CellType, height=7,width=7, dpi=300, filename=paste0("Fig2B_OutUMAP_MainCellType_Rat_snRNAseq_YesLabelSubset.pdf"), useDingbats=FALSE)

##### $$$$$ ###### $$$$$$ ###### Figure 2C. UMAP Feature plot by marker gene expression.   ##### $$$$$ ###### $$$$$ #####
# MyFeaturePlot <- Seurat::FeaturePlot(SeuratObject_Rat_snRNAseq,raster=TRUE, cols=c("lightgrey","darkblue"), order=TRUE, pt.size=1.8, 
#                                      features=c("Epcam","Krt18" )) #"KRT18","CD14", "PDGFRB","ACTA2", "CD3D","CD68","MS4A1", "PECAM1","IL7R",'CD8A'))   #  raster.dpi = c(512, 512)   
# ggsave(MyFeaturePlot, height=12.1,width=18, dpi=300, filename=paste0("Fig2C_OutFeatureUMAP_Rat_snRNAseq_ByMainMarkergeneExpression.pdf"), useDingbats=FALSE)
##### ----------- ============== $$$$$$$$$$$$$  ##### ----------- ============== $$$$$$$$$$$$$  ##### ----------- ============== $$$$$$$$$$$$$  

## Metadata has the cluster IDs already. The last column, "seurat_clusters" has the cluster IDS
SeuratObject_Rat_snRNAseq$CaseID <- SeuratObject_Rat_snRNAseq$orig.ident; table(SeuratObject_Rat_snRNAseq$CaseID )
table(SeuratObject_Rat_snRNAseq@active.ident)


MyMainMetadata <- SeuratObject_Rat_snRNAseq@meta.data %>% tibble::rownames_to_column("CellID_RatsnRNAseq") %>% 
  dplyr::select(-c(orig.ident,nCount_RNA,nFeature_RNA,percent.mt,SCT_snn_res.1.5));dim(MyMainMetadata) # 94129  8
table(MyMainMetadata$seurat_clusters); table(MyMainMetadata$AgeGroup) # Aged: 50792,  Young: 43337
table(MyMainMetadata$CellTypeByMarker_RatsnRNAseq)
MyMainMetadata[1:2,]
#                   CellID_RatsnRNAseq percent.rbp nCount_SCT nFeature_SCT AgeGroup seurat_clusters             CaseID CellTypeByMarker_RatsnRNAseq
# 1 Lee_021924_Nuclei1_AAACCCAAGAATTGTG-1           0       4323         2868     Aged               5 Lee_021924_Nuclei1             CancerEpithelial
# 2 Lee_021924_Nuclei1_AAACCCAAGACATAAC-1           0       3005         1945     Aged              15 Lee_021924_Nuclei1                   ImmuneCell

### Violin plot by marker genes. 
MyViolinplot_Main <- Seurat::VlnPlot(object=SeuratObject_Rat_snRNAseq, features=c("Epcam", "Krt18", "Mylk","Myh11", "Ptprc",  "Stat4",   "Cd74", "Lyz2",  "Col1a1","Pecam1"), stack=TRUE,flip=TRUE ) # c("CPA3", "PECAM1","EPCAM","KRT19","PDGFRB","ACTA2","CD68","MS4A1","CD3D")
ggsave(MyViolinplot_Main, height=6,width=6, dpi=300, filename=paste0("OutViolin_MainCellTypeMarker_HarmoneyRatsnRNAseq_Res1.0PC30KP15_Subtype.pdf"), useDingbats=FALSE)

### Store CellID, CaseID, CellTypeAnno data to txt file. 
# MedtaData <- SeuratObject_Rat_snRNAseq@meta.data[, c("CaseID","CellTypeByMarker_RatsnRNAseq")]; dim(MedtaData); MedtaData[1:2,] # 30959     3
MyMainMetadata_Proc <- MyMainMetadata %>% dplyr::select(c(CellID_RatsnRNAseq, CaseID,AgeGroup,seurat_clusters,CellTypeByMarker_RatsnRNAseq))
colnames(MyMainMetadata_Proc)[1] <- "CellID"; MyMainMetadata_Proc[1:2,]; dim(MyMainMetadata_Proc)
# MedtaData_CellID <- MyMainMetadata_Proc %>% tibble::rownames_to_column("CellID"); dim(MedtaData_CellID)
fwrite(MyMainMetadata_Proc, file="Metadata_RatsnRNAseq_YoungElderly_94129c.txt", row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t") # 94129     5
```


#### =========== ####  Section 5. Subset Myeloid cells and find macrophage cells by marker genes expression  #### =========== ####

```{r}

################################################################################
## Step5a. Subset Myeloid and define cell subtypes
################################################################################
SeuratObject_Rat_snRNAseq_Myeloid <- base::subset(x=SeuratObject_Rat_snRNAseq, idents =c("Myeloid")) ##    subset(x=Seurat_object, idents = c("Myeloid")) # Subset by active.ident
length(SeuratObject_Rat_snRNAseq_Myeloid@active.ident)  # 4822

HighVarGene <- SeuratObject_Rat_snRNAseq_Myeloid@assays$SCT@var.features;length(HighVarGene)
SeuratObject_Rat_snRNAseq_Myeloid <- Seurat::ScaleData(SeuratObject_Rat_snRNAseq_Myeloid, features=HighVarGene)
SeuratObject_Rat_snRNAseq_Myeloid<- Seurat::RunPCA(SeuratObject_Rat_snRNAseq_Myeloid, verbose = FALSE)  #

SeuratObject_Rat_snRNAseq_Myeloid <- Seurat::FindNeighbors(SeuratObject_Rat_snRNAseq_Myeloid, reduction = "pca", dims = 1:30)
SeuratObject_Rat_snRNAseq_Myeloid <- Seurat::FindClusters(SeuratObject_Rat_snRNAseq_Myeloid, resolution = 1.0)  # Res0.4: 12 clusters # Res0.8, 18 clusters # Res0.1  6 clusters
table(SeuratObject_Rat_snRNAseq_Myeloid@active.ident)

saveRDS(SeuratObject_Rat_snRNAseq_Myeloid, file="SeuratObj_Rat_snRNAseq_Myeloid.rds")

## By TSNE
SeuratObject_Rat_snRNAseq_Myeloid <- Seurat::RunTSNE(SeuratObject_Rat_snRNAseq_Myeloid, dims = 1:30)  # lower dims numbers will give less number of clusters image (less white space)
MyDimplot_ERpos_Myeloid_TSNE <- Seurat::DimPlot(SeuratObject_Rat_snRNAseq_Myeloid, pt.size=1, label.size=10,reduction="tsne", label=TRUE)
ggsave(MyDimplot_ERpos_Myeloid_TSNE, height=8,width=9, dpi=300, filename=paste0("OutTSNEplot_RatsnRNAseq_MacrophageDC_Res1.2_21clst_FromWholeRes1.0PC30KP15_ERp.pdf"), useDingbats=FALSE)

saveRDS(SeuratObject_Rat_snRNAseq_Myeloid, file="SeuratObj_Rat_snRNAseq_Myeloid_TSNE.rds")

################################################################################
## Step5b. Dotplot by Myeloid markers
################################################################################
# MyFeature_Myeloid <- c( "IL1B","CSF1R","S100A9","FCGR3A",   "CD14","SIGLEC1","CXCL10","EGR1","CD68","FCGR1A","FABP5","APOE",   "CD80","CD86","CD163","MRC1","MSR1",
#                         "CLEC10A","THBD","CD1C","ITGAX","HLA-DRB1","LAMP3","CLEC9A","GZMB","FLT3",     "IL3RA","SELL","IRF7","TYROBP", ### I used these 30 genes for cell type assignment
#                         "LILRA4","CXCR3","CLEC4C",   "IL10","CD40",   "CCR2","CCL2","CCL18","MMP9","CX3CR1","MT1G","SLC2A1","LYVE1","LYZ",   "ACE", "ADGRE1",  "TEK")
# # "ADGRE1","TEK", don't have expression.
# MyFeature_Myeloid_Cheng <- c("LST1","LILRB2","FCN1","CD14", "PTPRC","CD93","S100A8","S100A9","CD86",     "C1QA","C1QC","SEPP1","PLTP","EREG","NLRP3","CCL4","CD68","FCGR3A","MRC1",   "IDO1","CCR7","FSCN1","SELL","IRF7")
MyFeature_Myeloid <- c("S100a8","Adgre1","Ace","Csf1r",
                        "Il1b","S100a9","Fcgr3a",   "Cd14","Siglec1","Cxcl10","Egr1","Cd68","Fcgr1a","Fabp5","Apoe",   "Cd80","Cd86","Cd163","Mrc1","Msr1","Clec10a",   "Thbd","Itgax","Lamp3","Clec9a","Gzmb","Flt3", 
                        "Il3ra","Sell","Irf7")

MyDotPlot_Myeloid <- Seurat::DotPlot(SeuratObject_Rat_snRNAseq_Myeloid, features=MyFeature_Myeloid)+ theme(axis.text.x=element_text(vjust=0.6, angle=90)) +
  scale_colour_gradient2(low = "#1515FA", mid = "#FFFAE2", high = "#ff0000") + coord_flip()    # scale_color_viridis_c()    #ffe272
ggsave(MyDotPlot_Myeloid, height=8,width=8, dpi=300, filename=paste0("OutDotplot_ByMyeloidMarker_RatsnRNAseq_Res1.0_25clst_FromWholecellRes1.0.pdf"), useDingbats=FALSE)

## Metadata has the cluster IDs already. The last column, "seurat_clusters" has the cluster IDS
MyMetadata <- SeuratObject_Rat_snRNAseq_Myeloid@meta.data %>% tibble::rownames_to_column("CellID_") %>%
  dplyr::select(-c(orig.ident,nCount_RNA,nFeature_RNA,percent.mt,SCT_snn_res.1));dim(MyMetadata) # 4822   9
table(MyMetadata$seurat_clusters)
#   0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24 
# 542 538 448 359 338 243 235 232 232 211 197 167 152 152 135 113 101  81  79  73  58  42  33  32  29 

# ## Macrophage assignment by Wu et al. paper.
# MyMetadata_Mac <- MyMetadata %>% dplyr::filter(CellTypeMinor=="Macrophage", ); dim(MyMetadata_Mac)  # 1578  27
# table(MyMetadata_Mac$seurat_clusters)  # 0, 7, 8, 10, 16 are macrophage, maybe 9 is not macrophage

## from whole cell res 1.0  
new.cluster.ids_Myeloid <- c("Monocyte","Monocyte","Monocyte","Macrophage","Monocyte","M2_Macrophage",            "M2_Macrophage","M2_Macrophage","Monocyte","Macrophage","Monocyte",
                             "M2_Macrophage","Monocyte","Monocyte","Macrophage","Monocyte",                        "Monocyte","Monocyte","M1_Macrophage","Monocyte","Monocyte",
                             "Monocyte", "DendriticCell","Monocyte","Monocyte")  
names(new.cluster.ids_Myeloid) <- levels(SeuratObject_Rat_snRNAseq_Myeloid)
SeuratObject_Rat_snRNAseq_Myeloid <- Seurat::RenameIdents(SeuratObject_Rat_snRNAseq_Myeloid, new.cluster.ids_Myeloid)
table(SeuratObject_Rat_snRNAseq_Myeloid@active.ident)
# Monocyte    Macrophage M2_Macrophage M1_Macrophage      DendriticCell 
# 3128           705           877            79            33 

MeyloidSubsetMetadata <- SeuratObject_Rat_snRNAseq_Myeloid@meta.data[, c("CaseID","seurat_clusters" )]
MeyloidSubtypeDF <- data.frame(SeuratObject_Rat_snRNAseq_Myeloid@active.ident); colnames(MeyloidSubtypeDF)[1]<-"MyeloidSubtype"
MeyloidClusterNumbSubtype <- cbind(MeyloidSubsetMetadata,MeyloidSubtypeDF )
# fwrite(MeyloidClusterNumbSubtype, file="MeyloidClusterNumbSubtype_ERp.txt", row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)

## By TSNE
SeuratObject_Rat_snRNAseq_Myeloid <- Seurat::RunTSNE(SeuratObject_Rat_snRNAseq_Myeloid, dims = 1:30)  # lower dims numbers will give less number of clusters image (less white space)
MyDimplot_ERpos_Myeloid_Rename_TSNE <- Seurat::DimPlot(SeuratObject_Rat_snRNAseq_Myeloid, pt.size=1, label.size=8, reduction="tsne", label=TRUE)
# ggsave(MyDimplot_ERpos_Myeloid_Rename_TSNE, height=8,width=9, dpi=300, filename=paste0("OutTSNE_ByMyeloidMarker_RatsnRNAseq_Res1.2_21clst_FromWholeRes1.0PC30KP15_ERp.pdf"), useDingbats=FALSE)

# ## Cell type annotation by Wu et al. and Bassez et al.
# SeuratObject_Rat_snRNAseq_Myeloid$CellType_WuBassez <- ifelse(is.na(SeuratObject_Rat_snRNAseq_Myeloid$cellType),paste0(SeuratObject_Rat_snRNAseq_Myeloid$CellTypeMinor, "_Wu"), paste0(SeuratObject_Rat_snRNAseq_Myeloid$cellType, "_Bassez" ) ) #
# MyDimplot_ERpos_Myeloid_Rename_TSNE_Wu <- Seurat::DimPlot(SeuratObject_Rat_snRNAseq_Myeloid, pt.size=1, label.size=8, reduction="tsne", label=TRUE, group.by="CellType_WuBassez")
# #ggsave(MyDimplot_ERpos_Myeloid_Rename_TSNE_Wu, height=8,width=9, dpi=300, filename=paste0("OutTSNE_ByMyeloidMarker_RatsnRNAseq_WuBassezCellType.pdf"), useDingbats=FALSE)
# table(SeuratObject_Rat_snRNAseq_Myeloid@active.ident, SeuratObject_Rat_snRNAseq_Myeloid$CellType_WuBassez)

MyFeaturePlot_Myeloid <- Seurat::FeaturePlot(SeuratObject_Rat_snRNAseq_Myeloid,raster=FALSE,reduction="tsne", pt.size=0.5,features = c("Cd68","Fcgr3a","Cd14","Cd86","Mrc1","Csf1r","Sell","Irf7"))   #  raster.dpi = c(512, 512)
# ggsave(MyFeaturePlot_Myeloid, height=10,width=10, dpi=300, filename=paste0("OutFeatureUMAP_Myeloid_Res1.0PC30KP15_Res1.2_ERp.pdf"), useDingbats=FALSE)

##### $$$$$ ###### $$$$$$ ###### Supplementary Figure 1B. Feature violin plot by marker gene expression for Myeloid   ##### $$$$$ ###### $$$$$ #####
MyViolinplot <- Seurat::VlnPlot(object = SeuratObject_Rat_snRNAseq_Myeloid, features = c("Cd68","Fcgr3a","Cd14","Cd86",  "Mrc1","Msr1",  "Csf1r","Adgre1","Fcgr3a", "Gzmb", "Fabp5"),pt.size=2, stack=TRUE, flip=TRUE )
ggsave(MyViolinplot, height=5,width=5, dpi=300, filename=paste0("FigS1B_OutViolinplot_MyeloidMarkerExp.pdf"), useDingbats=FALSE)
##### ----------- ============== $$$$$$$$$$$$$  ##### ----------- ============== $$$$$$$$$$$$$  ##### ----------- ============== $$$$$$$$$$$$$
```

```{r}
################################################################################
## Step5c. Store the cell type identity as a data.frame
################################################################################
##### Store the cell type identity as a data.frame
CellTypeMyeloidDataFrame <- data.frame((SeuratObject_Rat_snRNAseq_Myeloid@active.ident)); dim(CellTypeMyeloidDataFrame) # 4822   1
colnames(CellTypeMyeloidDataFrame)[1]<-"MyeloidType";
CellTypeMyeloidDataFrame_Proc <- CellTypeMyeloidDataFrame %>% tibble::rownames_to_column("CellID_RatsnRNAseq"); head(CellTypeMyeloidDataFrame_Proc); table(CellTypeMyeloidDataFrame_Proc$MyeloidType)
# Monocyte    Macrophage M2_Macrophage M1_Macrophage            DendriticCell 
#     3128           705           877            79            33 
```

```{r}
################# ================ ################# ================ ################# ================ ################# ================ ################# ================
#### =========== ####   Section 6. Subset NKT cells and find cell sub types by marker genes expression #### =========== ####
################# ================ ################# ================ ################# ================ ################# ================ ################# ================
################################################################################
## Step6a. Subset NKT cells
################################################################################
SeuratObject_Rat_snRNAseq_NKT <- base::subset(x=SeuratObject_Rat_snRNAseq, idents =c("NKTcell")) # Subset by active.ident
length(SeuratObject_Rat_snRNAseq_NKT@active.ident)  # 3608

# # Normalize => Find Variable Features
HighVarGene <- SeuratObject_Rat_snRNAseq_NKT@assays$SCT@var.features;length(HighVarGene)  # 3000
SeuratObject_Rat_snRNAseq_NKT <- Seurat::ScaleData(SeuratObject_Rat_snRNAseq_NKT, features=HighVarGene)
SeuratObject_Rat_snRNAseq_NKT<- Seurat::RunPCA(SeuratObject_Rat_snRNAseq_NKT, verbose = FALSE)  #features = VariableFeatures(object=SeuratObject_Rat_snRNAseq_NKT), npcs=10)  # large npcs number takes long.   50 is better than 10

SeuratObject_Rat_snRNAseq_NKT <- Seurat::FindNeighbors(SeuratObject_Rat_snRNAseq_NKT, reduction="pca", dims = 1:30)
SeuratObject_Rat_snRNAseq_NKT <- Seurat::FindClusters(SeuratObject_Rat_snRNAseq_NKT, resolution=1.0)
table(SeuratObject_Rat_snRNAseq_NKT@active.ident)

saveRDS(SeuratObject_Rat_snRNAseq_NKT, file="SeuratObj_HarmonyRatsnRNAseqs_NKTcell_FromWholeRes1.0PC30KP15_Res1.0_ERp.rds")

SeuratObject_Rat_snRNAseq_NKT <- Seurat::RunTSNE(SeuratObject_Rat_snRNAseq_NKT, dims = 1:30)  # lower dims numbers will give less number of clusters image (less white space)
MyDimplot_ERpos_Tcell_TSNE <- Seurat::DimPlot(SeuratObject_Rat_snRNAseq_NKT, pt.size=1, label.size=8,reduction="tsne", label=TRUE)
# ggsave(MyDimplot_ERpos_Tcell_TSNE, height=8,width=8, dpi=300, filename=paste0("OutTSNEplot_RatsnRNAseq_Tcell_Res1.0_31clst_FromWholeRes1.0PC30KP15_ERp.pdf"), useDingbats=FALSE)

#################################################################################
## Step6b. Dotplot by NK, CD4+T, CD8+T cell markers
#################################################################################
# Marker_Tcell <- c("FCGR3A","IL2RB","KLRB1","KLRC1","KLRD1","KLRF1","KLRK1","NCR3","NCR1","FGFBP2","AREG","XCL1","KIR2DL4","NCAM1",
#                  "CD8A","CD8B","GZMA","GZMB",  "LAG3","PDCD1","CTLA4",   "TNFRSF4","BTLA","CD40LG","STAT4","STAT1","FOXP3","CD4","CCR7",    "IL7R","IL2RA","PTPRC","CD3D","CD27")   # ,

Marker_Tcell <- c("Fcgr3a","Il2rb","Klrb1","Klrc1","Klrd1","Klrk1","Ncr3","Ncr1","Areg","Xcl1","Ncam1",
                  "Cd8a","Cd8b","Gzma","Gzmb",  "Lag3","Pdcd1","Ctla4",   "Tnfrsf4","Btla","Cd40lg","Stat4","Foxp3","ENSRNOG00000071219","Ccr7",    "Il7r","Il2ra","Ptprc","Cd3d","Cd27")   # ENSRNOG00000071219 is Cd4

MyDotPlot_Tcell <- Seurat::DotPlot(SeuratObject_Rat_snRNAseq_NKT, features=Marker_Tcell)+ theme(axis.text.x=element_text(vjust=0.6, angle=90)) +
  scale_colour_gradient2(low = "#1515FA", mid = "#FFFAE2", high = "#ff0000") + coord_flip()    # scale_color_viridis_c()    #ffe272
# ggsave(MyDotPlot_Tcell, height=8,width=8, dpi=300, filename=paste0("OutDotplot_ByTcellMarker_RatsnRNAseq_Res0.8_17clst_FromWholecellRes1.5_ERp.pdf"), useDingbats=FALSE)

## Metadata has the cluster IDs already. The last column, "seurat_clusters" has the cluster IDS
MyMetadata <- SeuratObject_Rat_snRNAseq_NKT@meta.data %>% tibble::rownames_to_column("CellID") %>%
  dplyr::select(-c(orig.ident,nCount_RNA,nFeature_RNA,percent.mt,SCT_snn_res.1));dim(MyMetadata) # 16575    14
table(MyMetadata$seurat_clusters)

#### =========== New cell type assignment   ============= ###########
new.cluster.ids_Tcell <- c("NKcell","CD8Tcell","CD4Tcell","CD8Tcell","CD4Tcell","NKcell",          "NaiveTcell","CD4Tcell","NKcell","Treg","CD8Tcell",
                           "CD8Tcell","CD8Tcell","NKcell","NKcell","NKcell",                    "NaiveTcell","CD8Tcell","NaiveTcell")
names(new.cluster.ids_Tcell) <- levels(SeuratObject_Rat_snRNAseq_NKT)
SeuratObject_Rat_snRNAseq_NKT <- Seurat::RenameIdents(SeuratObject_Rat_snRNAseq_NKT, new.cluster.ids_Tcell)
table(SeuratObject_Rat_snRNAseq_NKT@active.ident); sum(table(SeuratObject_Rat_snRNAseq_NKT@active.ident))  # 16956
# NKcell   CD8Tcell   CD4Tcell NaiveTcell       Treg 
# 1126       1183        827        306        166 
# saveRDS(SeuratObject_Rat_snRNAseq_NKT, "SeuratObject_Rat_snRNAseq_NKTSubtyping.rds")

## subset to exclude 'Remove' type cells  =========== Just for this TSNE plot, let's exclude 'NotAvail' cells
# SeuratObject_Rat_snRNAseq_NKT <- base::subset(x=SeuratObject_Rat_snRNAseq_NKT, idents = c("NotAvail"), invert=TRUE)  # invert=TRUE willl exclude the ident cells.
# table(SeuratObject_Rat_snRNAseq_NKT@active.ident)

SeuratObject_Rat_snRNAseq_NKT <- Seurat::RunTSNE(SeuratObject_Rat_snRNAseq_NKT, dims = 1:30)  # lower dims numbers will give less number of clusters image (less white space)
MyDimplot_ERpos_Tcell_Rename_TSNE <- Seurat::DimPlot(SeuratObject_Rat_snRNAseq_NKT,label.size=8, pt.size=1, reduction="tsne", label=FALSE)
# ggsave(MyDimplot_ERpos_Tcell_Rename_TSNE, height=8,width=8.5, dpi=300, filename=paste0("OutTSNE_ByTcellMarker_RatsnRNAseq_Res1.0_31clst_FromWholeRes1.5PC30KP15_ERp.pdf"), useDingbats=FALSE)

##### $$$$$ ###### $$$$$$ ###### Supplementary Figure 1C. Feature violin plot by marker gene expression for NK/T cell   ##### $$$$$ ###### $$$$$ #####
SeuratObject_Rat_snRNAseq_NKT@active.ident <- factor(SeuratObject_Rat_snRNAseq_NKT@active.ident, levels=c("NKcell","NaiveTcell","CD4Tcell","CD8Tcell","Treg"))
MyViolinplot_NKTcell <- Seurat::VlnPlot(object = SeuratObject_Rat_snRNAseq_NKT, features = c("Klrd1","Areg","Gzmb","Gzma", "Xcl1", "Il7r","Il2ra","ENSRNOG00000071219","Cd8a","Cd8b","Foxp3","Cd27"), stack=TRUE, flip=TRUE) # c("AREG","KLRB1", "XCL1", "IL7R","IL2RB","CD4","CD8A","CD8B","FOXP3")
ggsave(MyViolinplot_NKTcell, height=4.5,width=8, dpi=300, filename=paste0("FigS1C_OutViolin_NKTcellMarkerExp.pdf"), useDingbats=FALSE)
##### ----------- ============== $$$$$$$$$$$$$  ##### ----------- ============== $$$$$$$$$$$$$  ##### ----------- ============== $$$$$$$$$$$$$

```

```{r}
#################################################################################
## Step6c. Store the cell type identity as a data.frame
#################################################################################
CellTypeTcellDataFrame <- data.frame((SeuratObject_Rat_snRNAseq_NKT@active.ident)); dim(CellTypeTcellDataFrame) # 3608   1
colnames(CellTypeTcellDataFrame)[1]<-"TcellType";
CellTypeTcellDataFrame_Proc <- CellTypeTcellDataFrame %>% tibble::rownames_to_column("CellID_RatsnRNAseq"); head(CellTypeTcellDataFrame_Proc)
table(CellTypeTcellDataFrame_Proc$TcellType); dim(CellTypeTcellDataFrame_Proc)  # 3608
# NKcell   CD8Tcell   CD4Tcell NaiveTcell       Treg 
# 1126       1183        827        306        166 

```

# Section 7.  Annotate Myeloid cell subtypes and NKT cell subtypes, and calculate Macrophage fraction per samples.    #### =========== ####

## Step7a. Annotate Myeloid cell subtypes and NKT cell subtypes

```{r}
CellTypeClusterDataFrame_Proc <- MyMainMetadata[, c("CellID_RatsnRNAseq","CellTypeByMarker_RatsnRNAseq")]
CellTypeClusterMyeloidTcell <- dplyr::left_join(CellTypeClusterDataFrame_Proc, CellTypeMyeloidDataFrame_Proc) %>% dplyr::left_join(CellTypeTcellDataFrame_Proc)
dim(CellTypeClusterMyeloidTcell) # [1] 94129     3
CellTypeClusterMyeloidTcell[1:2,]
table(CellTypeClusterMyeloidTcell$CellTypeByMarker_RatsnRNAseq)
# CancerEpithelial    Myoepithelial          Myeloid          NKTcell       Fibroblast    DendriticCell      Endothelial 
#         77318             5550             4822             3608             1230              838              763 

CellTypeClusterMyeloidTcell$ByMainMarker_MyeloidTcell <- ifelse(CellTypeClusterMyeloidTcell$CellTypeByMarker_RatsnRNAseq == "Myeloid",
                                                                as.vector(CellTypeClusterMyeloidTcell$MyeloidType), ifelse(CellTypeClusterMyeloidTcell$CellTypeByMarker_RatsnRNAseq == "NKTcell",
                                                                as.vector(CellTypeClusterMyeloidTcell$TcellType), as.vector(CellTypeClusterMyeloidTcell$CellTypeByMarker_RatsnRNAseq)))
table(CellTypeClusterMyeloidTcell$ByMainMarker_MyeloidTcell); sum(table(CellTypeClusterMyeloidTcell$ByMainMarker_MyeloidTcell))  # 94129
# CancerEpithelial         CD4Tcell         CD8Tcell    DendriticCell      Endothelial       Fibroblast    M1_Macrophage    M2_Macrophage       Macrophage         Monocyte    Myoepithelial       NaiveTcell           NKcell             Treg 
# 77318              827             1183              871              763             1230               79              877              705             3128             5550              306             1126              166 
################################################################################
## Step7b. Add new metadata of new cell subtypes to Seurat object
################################################################################
SeuratObject_Rat_snRNAseq <- Seurat::AddMetaData(object=SeuratObject_Rat_snRNAseq, metadata=c(CellTypeClusterMyeloidTcell$ByMainMarker_MyeloidTcell), col.name=c("CellTypeMacroTcell_RatsnRNAseq"))
table(SeuratObject_Rat_snRNAseq$CellTypeMacroTcell_RatsnRNAseq) ## New cell type on 20230612
# CancerEpithelial         CD4Tcell         CD8Tcell    DendriticCell      Endothelial       Fibroblast    M1_Macrophage    M2_Macrophage       Macrophage         Monocyte    Myoepithelial       NaiveTcell           NKcell             Treg 
#         77318              827             1183              871              763             1230               79              877              705             3128             5550              306             1126              166 
SeuratObject_Rat_snRNAseq <- Seurat::SetIdent(SeuratObject_Rat_snRNAseq, value=SeuratObject_Rat_snRNAseq$CellTypeMacroTcell_RatsnRNAseq)
table(SeuratObject_Rat_snRNAseq@active.ident)

saveRDS(SeuratObject_Rat_snRNAseq, file="SeuratObj_Harmony_Rat_snRNAseq_CelltypeDefined.rds")

################################################################################
## Step7c. Remove "NotAvail" cells
################################################################################
# SeuratObject_Rat_snRNAseq <- base::subset(x=SeuratObject_Rat_snRNAseq, idents = c("NotAvail"), invert=TRUE)  # invert=TRUE willl exclude the ident cells.
# table(SeuratObject_Rat_snRNAseq@active.ident)
table(SeuratObject_Rat_snRNAseq$CellTypeMacroTcell_RatsnRNAseq); sum(table(SeuratObject_Rat_snRNAseq$CellTypeMacroTcell_RatsnRNAseq)) # 94129
# CancerEpithelial         CD4Tcell         CD8Tcell    DendriticCell      Endothelial       Fibroblast    M1_Macrophage    M2_Macrophage       Macrophage         Monocyte    Myoepithelial       NaiveTcell           NKcell             Treg 
# 77318              827             1183              871              763             1230               79              877              705             3128             5550              306             1126              166 

CellTypeClusterMyeloidTcell_Rmv <- CellTypeClusterMyeloidTcell %>% dplyr::filter(!ByMainMarker_MyeloidTcell %in% c("NotAvail")); dim(CellTypeClusterMyeloidTcell_Rmv) # 94129 5 

table(CellTypeClusterMyeloidTcell_Rmv$ByMainMarker_MyeloidTcell)  # This has NO "NotAvail)
saveRDS(SeuratObject_Rat_snRNAseq, file="SeuratObj_Harmony_Rat_snRNAseq_CelltypeDefined_RmvNotAvail.rds")

################# ======================= MetaData after cell type assignment.  ======================= #################
MyMainMetadata_AfterCellType <- SeuratObject_Rat_snRNAseq@meta.data %>% data.frame; dim(MyMainMetadata_AfterCellType) # 94129    13
MyMainMetadata_AfterCellType_Proc <- MyMainMetadata_AfterCellType %>% tibble::rownames_to_column("CellID") %>% dplyr::select(CaseID,CellID, AgeGroup, seurat_clusters,CellTypeByMarker_RatsnRNAseq, CellTypeMacroTcell_RatsnRNAseq)
saveRDS(MyMainMetadata_AfterCellType_Proc, "MyMainMetadata_AfterCellType_Proc.rds")
################# ======================= ################# ======================= ################# ======================= 


### Violin plot by marker genes per main cell types
SeuratObject_Rat_snRNAseq$CellTypeByMarker_RatsnRNAseq <- factor(SeuratObject_Rat_snRNAseq$CellTypeByMarker_RatsnRNAseq, levels=c("CancerEpithelial","Myoepithelial","Endothelial","Fibroblast","Myeloid","DendriticCell","NKTcell"))
                                                                                                                                      
MyViolinplot_MainType <- Seurat::VlnPlot(object=SeuratObject_Rat_snRNAseq, features=c("Epcam", "Krt18","Acta2", "Mylk","Myh11", "Ptprc",  "Stat4",   "Cd74", "Lyz2","Csf1r",  "Fabp5","Msr1","Mrc1", "Cd86","Cd68",  "Pdgfrb", "Col1a1","Pecam1"),
                                       group.by="CellTypeByMarker_RatsnRNAseq", stack=TRUE,flip=TRUE ) #  group.by="CellTypeMacroTcell_RatsnRNAseq"
ggsave(MyViolinplot_MainType, height=6,width=7, dpi=300, filename=paste0("OutViolin_MainCellTypeMarker_HarmoneyRatsnRNAseq_Res1.0PC30KP15_20240626.pdf"), useDingbats=FALSE)

### Violin plot by marker genes per sub celltypes   <<<== I should load "SeuratObject_Rat_snRNAseq" first at line 491
SeuratObject_Rat_snRNAseq$CellTypeMacroTcell_RatsnRNAseq <- factor(SeuratObject_Rat_snRNAseq$CellTypeMacroTcell_RatsnRNAseq, levels=c("CancerEpithelial","Myoepithelial","Endothelial","Fibroblast","Monocyte",
                                                                                                                                      "Macrophage","M1_Macrophage","M2_Macrophage","DendriticCell","NKcell","CD4Tcell","CD8Tcell","Treg","NaiveTcell"))
MyViolinplot_MainSub <- Seurat::VlnPlot(object=SeuratObject_Rat_snRNAseq, features=c("Epcam", "Krt18","Acta2", "Mylk","Myh11", "Ptprc",  "Stat4",   "Cd74", "Lyz2","Csf1r",  "Fabp5","Msr1","Mrc1", "Cd86","Cd68",  "Pdgfrb", "Col1a1","Pecam1"),
                                        group.by="CellTypeMacroTcell_RatsnRNAseq", stack=TRUE,flip=TRUE ) #  group.by="CellTypeMacroTcell_RatsnRNAseq" or "CellTypeByMarker_RatsnRNAseq"
ggsave(MyViolinplot_MainSub, height=6.5,width=6.5, dpi=300, filename=paste0("OutViolin_SubCellTypeMarker_HarmoneyRatsnRNAseq_Res1.0PC30KP15_20240624.pdf"), useDingbats=FALSE)

### Dotplot visualizing averaged expression of cannonical markers in cell clusters. .
# MyMainMarker <- c("PECAM1","RAMP2","FLT1","CLDN5",   "EPCAM","KRT19","KRT18","CD24",  "PDGFRB","C1R","DCN","COL1A1",   "ACTA2",  "TPSB2","TPSAB1","CPA3",
#                   "CD68","LYZ","TYROBP",   "LILRA4","LAMP3",  "CLEC9A","STAT1","IRF7","CD83","MS4A1","MZB1","CD79A",  "CD8A","CD2","CD3E","CD3D","CD3G","IL7R") # I should include "CD79A","CD19"  B cell marker.
# MyDotPlot_MainType <- Seurat::DotPlot(SeuratObject_Rat_snRNAseq, features=MyMainMarker, group.by="CellTypeMacroTcell_RatsnRNAseq")+ theme(axis.text.x=element_text(vjust=0.6, size=15,angle=0), axis.text.y=element_text(vjust=0.6, size=15,angle=0)) +
#   scale_colour_gradient2(low = "#1515FA", mid = "#FFFAE2", high = "#ff0000") + coord_flip()    # scale_color_viridis_c()    #ffe272
# # ggsave(MyDotPlot_MainType, height=10,width=15, dpi=300, filename=paste0("OutDotplot_ByMainMarker_HarmoneyRatsnRNAseq_Res1.0_40Cluster_MoreMarker_20240222.pdf"), useDingbats=FALSE)

################ ================ ################# ================ ################# ================ ################# ================ ################# ================
#### =========== ####   Section 8.  Calculate cell type fraction per samples. #### =========== ####
###               Make dotplot of Macrophage fraction per sample. Make fraction barplot for all samples 
################# ================ ################# ================ ################# ================ ################# ================ ################# ================ 
################################################################################
## Step8a. Cell type fraction per samples. Make dotplot for macrophage fraction. 
## number of cells per cluster: https://github.com/satijalab/seurat/issues/2825
################################################################################
SeuratObject_Rat_snRNAseq <- Seurat::AddMetaData(object=SeuratObject_Rat_snRNAseq, metadata=c(SeuratObject_Rat_snRNAseq$CellTypeMacroTcell_RatsnRNAseq), col.name=c("CellTypeMacroTcell_RatsnRNAseq"))
identical(SeuratObject_Rat_snRNAseq$CellTypeByMarker_RatsnRNAseq, SeuratObject_Rat_snRNAseq$CellTypeMacroTcell_RatsnRNAseq) # 

CellTypeFraction<- with(SeuratObject_Rat_snRNAseq@meta.data, table(orig.ident, CellTypeMacroTcell_RatsnRNAseq)); CellTypeFraction[1:3,]; dim(CellTypeFraction)   # class is table.  #6  4
#             CellTypeMacroTcell_RatsnRNAseq
#         orig.ident           CancerEpithelial CD4Tcell CD8Tcell DendriticCell Endothelial Fibroblast M1_Macrophage M2_Macrophage Macrophage Monocyte Myoepithelial NaiveTcell NKcell  Treg
# Lee_021924_Nuclei1            13966       13        1            74         184        685             0           391        194      106          1283         17    408    26

####   To make input data for cell type fraction correlation heatmap  ##################
CellTypeFraction_DF <- CellTypeFraction %>% matrix(ncol=ncol(CellTypeFraction)) %>% data.frame
colnames(CellTypeFraction_DF) <- colnames(CellTypeFraction)
rownames(CellTypeFraction_DF) <- rownames(CellTypeFraction)
head(CellTypeFraction_DF)

ProportionCalculation <- CellTypeFraction_DF/rowSums(CellTypeFraction_DF);  dim(ProportionCalculation) #  14  10
saveRDS(ProportionCalculation, "ProportionCalculation.rds")

library(ie2misc) # CRAN, for madstat()    ##   error: X11 library is missing: install XQuartz from www.xquartz.org
# cell_type <-"ImmuneCell" #  "Monocyte"
cell_type <- "Macrophage"
MacroFractMean <- mean(ProportionCalculation[,cell_type]*100); print(MacroFractMean)   # 9.426557
MacroFractMedian <- median(ProportionCalculation[,cell_type]*100);  print(MacroFractMedian)  # 6.619202
MacroFractSTD <- sqrt(var(ProportionCalculation[,cell_type]*100));  print(MacroFractSTD)  # 6.28527
MacroFractMAD <- sqrt(madstat(ProportionCalculation[,cell_type]*100));  print(MacroFractMAD) # 2.301107

dplyr::ntile(sort(ProportionCalculation[,cell_type])*100, 4) # 1 1 2 2 3 4
sort(ProportionCalculation[,cell_type])

##### ----------- ============== $$$$$$$$$$$$$  ##### ----------- ============== $$$$$$$$$$$$$  ##### ----------- ============== $$$$$$$$$$$$$  
##### $$$$$ ###### $$$$$$ ###### Supplementary Figure 1D. Dotplot to represent Macrophage infiltration in individual samples   ##### $$$$$ ###### $$$$$ #####
pdf(file="FigS1D_OutDotplot_MacrophageInfiltration_IndividualSample.pdf", width=5, height=5)
      plot(sort(ProportionCalculation[,cell_type])*100, pch=16, cex=1.3)
dev.off()
##### ----------- ============== $$$$$$$$$$$$$  ##### ----------- ============== $$$$$$$$$$$$$  ##### ----------- ============== $$$$$$$$$$$$$  

#################################################################################################
### Step8b. Make fraction barplot that represent each cell type fraction in different samples
#################################################################################################
NewIndex<-c()
#for (EachCase in names(sort(ProportionCalculation[,cell_type]))) { 
for (EachCase in rownames(ProportionCalculation)[order(ProportionCalculation[,cell_type])] ) { 
  Index<-print(grep(paste0(EachCase, "$"), rownames(ProportionCalculation) ) )
  NewIndex<-c(NewIndex, Index)
}

table(SeuratObject_Rat_snRNAseq$CellTypeMacroTcell_RatsnRNAseq)
# CancerEpithelial    Myoepithelial       ImmuneCell      StromalCell 
#       77318             5550             9268             1993          
# CancerEpithelial    Myoepithelial          Myeloid          NKTcell       Fibroblast    DendriticCell      Endothelial 
#       77318             5550             4822             3608             1230              838              763 

##### ----------- ============== $$$$$$$$$$$$$  ##### ----------- ============== $$$$$$$$$$$$$  ##### ----------- ============== $$$$$$$$$$$$$  
##### $$$$$ ###### $$$$$$ ######  Figure 1D. Barplot to represent cell type infiltration in individual samples   ##### $$$$$ ###### $$$$$ #####
FractionBarplot <- dittoBarPlot(object = SeuratObject_Rat_snRNAseq, var="CellTypeMacroTcell_RatsnRNAseq", group.by="orig.ident", x.reorder=NewIndex, var.labels.reorder=c(5,2,7,1,3,6,4),
                                color.panel=c("green","lightblue", "orange", "lightgrey","blue","red","yellow"))  #  x.reorder = seq(1,20,1)) x.reorder=c(5,3,4,1,2)
ggsave(FractionBarplot, height=5,width=8, dpi=300, filename=paste0("Fig2D_OutFractionBarplot_CellTypeInfiltration_Subtype.pdf"), useDingbats=FALSE)
##### ----------- ============== $$$$$$$$$$$$$  ##### ----------- ============== $$$$$$$$$$$$$  ##### ----------- ============== $$$$$$$$$$$$$  
```

```{r}
#################################################################################################
### Step8c. Make a table to record cell type fraction in different samples. 
#################################################################################################

# ProportionCalculationProcess <- ProportionCalculation %>% tibble::rownames_to_column("CaseID")
# CellTypeFraction_Tidyr <- ProportionCalculationProcess %>% tidyr::gather( key="CellType", value="CellFraction", colnames(ProportionCalculationProcess)[2]:colnames(ProportionCalculationProcess)[ncol(ProportionCalculationProcess)]); head(CellTypeFraction_Tidyr)  # use column names for tidyr range.
# table(CellTypeFraction_Tidyr$CaseID)

# ProportionCalculationProcess <- ProportionCalculation*100
# ProportionCalculationProcess$MacroPoorMidRich <- ifelse(ProportionCalculationProcess[, "ImmuneCell"] >= 7, "ImmuneRich", ifelse(ProportionCalculationProcess[, "ImmuneCell"] <= 2, "ImmunePoor","ImmuneMid"))
# ProportionCalculationMacroRatio <- ProportionCalculationProcess %>% tibble::rownames_to_column("CaseID"); colnames(ProportionCalculationMacroRatio)[ncol(ProportionCalculationMacroRatio)] <- "ImmuneCellRatio"
# table(ProportionCalculationMacroRatio$ImmuneCellRatio)
# ProportionCalculationMacroRatio$CaseID[ProportionCalculationMacroRatio$ImmuneCellRatio=="ImmuneMid"]   # CaseIDs of MacroMid 
# data.table::fwrite(ProportionCalculationMacroRatio, file="MacrophageInfiltrarion_IndividualSample.txt", col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)   ## Use this for celltype fraction correlation heatmap.

#################################################################################################
### Step8d. Cell type correlation pyramid plot.
#################################################################################################
# library(ComplexHeatmap)   #BiocManager::install("ComplexHeatmap")
# library(circlize)  # cran  for colorRamp2 function
# library(corrplot)
# options(rgl.useNULL = TRUE)  # https://stackoverflow.com/a/66127391/2554330
# 
# pathway.corplot_SH<-function(corMatrix,fontsize=0.8){
#   rownames(corMatrix)=gsub("_"," ",rownames(corMatrix))
#   colnames(corMatrix)=gsub("_"," ",colnames(corMatrix))
#   mycolor=colorRampPalette(c("darkblue","blue", "white","orange", "red","red"))(n=100)
#   corrplot(corMatrix, method="circle",type="lower",tl.col="black",col=mycolor,tl.cex=fontsize,tl.srt=30)
# }
# 
# CellTypeProportionData <- ProportionCalculationMacroRatio %>% tibble::column_to_rownames("CaseID") %>% dplyr::select(-c(MacrophageRatio ))
# 
# ## Calculate Spearman correlation.   === Calculate correlation with continuous values, not count data. 
# CellTypeProportionData_Correlation <- CellTypeProportionData [,c(1,  6,7,8,9,10,   2,3,4,5,11,12,13,14)]
# colnames(CellTypeProportionData_Correlation)
# 
# CellTypeProportionData_Cor <- cor(CellTypeProportionData_Correlation, method="spearman"); dim(CellTypeProportionData_Cor); CellTypeProportionData_Cor[1:2,] 
# 
# pdf(file="Fig2E_PyramidCorrelation_ByCellTypeFraction.pdf",width=6, height=6)
# pathway.plot=pathway.corplot_SH(corMatrix=CellTypeProportionData_Cor)
# graphics.off()
```



#  Section 9.  Scaled by DefaultAssay RNA - this is used to calculate correlation between gene expression and Macrophage fraction, or gene exp violin plot #### =========== ####
```{r}
#################################################################################################
## Step 9a. Scaled by DefaultAssay RNA 
#################################################################################################
library(patchwork) # plot_annotation function
library(gridExtra)  # for 'grid.arrange' function
library(ggrepel) # for geom_text_repel function. 

### Data scaling.
DefaultAssay(SeuratObject_Rat_snRNAseq) <- "RNA"
TotalGeneNumb<-nrow(SeuratObject_Rat_snRNAseq[["RNA"]]@features); print(TotalGeneNumb)  # 23098
SeuratObject_Rat_snRNAseq_Normal <- SeuratObject_Rat_snRNAseq %>% NormalizeData %>% FindVariableFeatures(selection.method="vst", nfeatures=TotalGeneNumb)  # nfeatures 20,000 is bettern than 5000
GeneSymb <- rownames(SeuratObject_Rat_snRNAseq_Normal@assays$RNA@features)
SeuratObject_Rat_snRNAseq_Scaled <- ScaleData(SeuratObject_Rat_snRNAseq_Normal, features=GeneSymb)  # This step requires a lot of memory. It may not run in your local computer. You will need a supercomputer.

saveRDS(SeuratObject_Rat_snRNAseq_Scaled, file="SeuratObject_Rat_snRNAseq_RNAScaled_Subset.rds")    
## SeuratObject_Rat_snRNAseq_Scaled <- readRDS("SeuratObject_Rat_snRNAseq_RNAScaled.rds")

#################################################################################################
## Step 9b. Gene expression violin plot. 
#################################################################################################
## subset by cell type
MyCellType <- names(table(SeuratObject_Rat_snRNAseq_Scaled$CellTypeByMarker_RatsnRNAseq)); print(MyCellType)
# [1] "CancerEpithelial" "Myoepithelial"    "ImmuneCell"       "StromalCell"    
# [1] "CancerEpithelial" "Myoepithelial"    "Myeloid"          "NKTcell"          "Fibroblast"       "DendriticCell"    "Endothelial"    

## Make violin plot per AgeGroup 
for(EachCellType in MyCellType) { # 
  # EachCellType <- MyCellType[3]; print(EachCellType)
  SeuratObj_ERpos_SubsetMyeloid <- subset(x=SeuratObject_Rat_snRNAseq_Scaled, idents=c(EachCellType)); table(SeuratObj_ERpos_SubsetMyeloid@active.ident) ## ImmuneCell: 9268
  SeuratObj_ERpos_SubsetMyeloid$AgeGroup <- factor(x=SeuratObj_ERpos_SubsetMyeloid$AgeGroup, levels=c("Young","Aged"))
  # SeuratObj_ERpos_SubsetMyeloid$orig.ident <- factor(x=SeuratObj_ERpos_SubsetMyeloid$orig.ident, levels=CaseID_SortByMacroFraction)
  
  MyFeature <- c("Cyp19a1","Esr1","Greb1","Hsd17b2","Hsd17b7","Pak4","Pgr") # Saa1 don't exist. 
  ##  ENSRNOG00000001164 is Tff1 but it doesn't exist. 
  
  MyViolinplot <- VlnPlot(object = SeuratObj_ERpos_SubsetMyeloid, features = c(MyFeature), group.by="AgeGroup", slot = 'scale.data')  # raster=TRUE requires ggrastr r packge.  # group.by = "orig.ident"
  OutFile <- paste0("OutVlnplot_EstrogenConversionGene_", EachCellType,  "_Subset_ByAgeGroup.pdf"); print(OutFile)
  ggsave(MyViolinplot, height=11,width=12, dpi=300, filename=OutFile, useDingbats=FALSE)  #  height=26,width=19,
}


## Sort CaseID by age
# SeuratObject_Rat_snRNAseq_Scaled$CaseID <- factor(SeuratObject_Rat_snRNAseq_Scaled$CaseID, levels=c(Clinicaldata_ProportionCal_Sort$CaseID))

## Make violin plot per each patients sorted by age. 
for(EachCellType in MyCellType) { # MyCellType_NoNASelc
  # EachCellType <- MyCellType[3]; print(EachCellType)
  SeuratObj_ERpos_SubsetMyeloid <- subset(x=SeuratObject_Rat_snRNAseq_Scaled, idents=c(EachCellType)); table(SeuratObj_ERpos_SubsetMyeloid@active.ident)
  #SeuratObj_ERpos_SubsetMyeloid$orig.ident <- factor(x=SeuratObj_ERpos_SubsetMyeloid$orig.ident, levels=CaseID_SortByMacroFraction)
  
  MyFeature <- c("Cyp19a1","Esr1","Greb1","Hsd17b2","Hsd17b7","Pak4","Pgr") # Saa1  don't exist.  

  
  MyViolinplot <- VlnPlot(object = SeuratObj_ERpos_SubsetMyeloid, features = c(MyFeature), group.by="CaseID", slot = 'scale.data')  # raster=TRUE requires ggrastr r packge.  # group.by = "orig.ident"
  OutFile <- paste0("OutVlnplot_EstrogenConversionGene_", EachCellType,  "_Subset_ByCaseID.pdf"); print(OutFile)
  ggsave(MyViolinplot, height=11,width=12, dpi=300, filename=OutFile, useDingbats=FALSE)  #  height=26,width=19,
}

```

This is to make CellPhoneDB or GSVA input files. 

```{r}
#################################################################################################
### =========== Step6. Make CellPhoneDB or GSVA input files   -  Baseline count data and cell metadata for YoungElderly  Rich vs Poor
#################################################################################################
## 3) Whole subjects   -- Can be used for GSVA and SCENIC too. 
# Remove zero count genes and save the count data for CellPhoneDB input
CountData_YoungElderlyWhole_CellType <- as.data.frame(SeuratObject_Rat_snRNAseq[["SCT"]]@counts); dim(CountData_YoungElderlyWhole_CellType)  # 18111 94129 cells
CountData_YoungElderlyWhole_CellType_NoLowCountGene <- CountData_YoungElderlyWhole_CellType[rowSums(CountData_YoungElderlyWhole_CellType)>10, ]; dim(CountData_YoungElderlyWhole_CellType_NoLowCountGene); #   150: 15979 30959  # 10: 17096 94129
## remove genes of non-proteincoding
CountData_YoungElderlyWhole_CellType_NoNonProtCoding <- CountData_YoungElderlyWhole_CellType_NoLowCountGene %>% dplyr::filter(!grepl("\\.", rownames(CountData_YoungElderlyWhole_CellType_NoLowCountGene)));
dim(CountData_YoungElderlyWhole_CellType_NoNonProtCoding) # 17096 94129
CountData_YoungElderlyWhole_CellType_NoLowCountGene <- CountData_YoungElderlyWhole_CellType_NoNonProtCoding[, colSums(CountData_YoungElderlyWhole_CellType_NoNonProtCoding)>1000 ]; dim(CountData_YoungElderlyWhole_CellType_NoLowCountGene); 
# 10/1000:  17096g 94122c
colnames(CountData_YoungElderlyWhole_CellType_NoLowCountGene) <- gsub("-","_",colnames(CountData_YoungElderlyWhole_CellType_NoLowCountGene))
fwrite(data.frame(CountData_YoungElderlyWhole_CellType_NoLowCountGene), "CountDataForGSVA_RatsnRNAseq_YoungElderlyWhole_MainCellType_17096g94122c.txt", col.names=TRUE, row.names=TRUE, sep="\t", quote=FALSE)  # 200/2000 19221g86138c
saveRDS(data.frame(CountData_YoungElderlyWhole_CellType_NoLowCountGene), file="CountDataForGSVA_RatsnRNAseq_YoungElderlyWhole_MainCellType_17096g94122c.rds")
```







