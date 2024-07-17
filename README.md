---
title: "Rat snRNA-seq data analysis"
format: html
editor: visual
R version: 4.4.1 (2024-16-14) -- "Race for Your Life"
---

### Rat snRNA-seq data analysis 

Systemic and breast chronic inflammation and hormone disposition promote a tumor-permissive locale for breast cancer in older women

![image](https://github.com/user-attachments/assets/a5a4ad70-20e1-4c36-976d-5cfe9c477cc4)

### Necessary input data files. (See Step2 below)

-   snRNA-seq processed count data 
-   Age group annotation data

## Preparation 1. Install and load the necessary R packages.
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

## Preparation 2. Set working directory and direct your input data files
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






