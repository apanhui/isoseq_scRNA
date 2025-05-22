## Pre-install software and R packages

### Software for Isoseq

- smrtlink_9.0.0.92188
- samtools-1.9
- cd-hit-v4.6.7
- diamond-0.9.25

### Software for Seurat and Cellchat

- R-4.1.3

### R packages for Seurat and Cellchat

- DoubletFinder-2.0.3
- Seurat-4.1.0
- CellChat-1.6.1
- harmony-0.1.0
- dplyr-1.0.9
- ggplot2-3.3.5
- randomcoloR-1.1.0.1
- patchwork-1.1.1
- future-1.30.0
- Cairo-1.6.0
- extrafont-0.19
- extrafontdb-1.0
- ggcor-0.9.8.1
- Rttf2pt1-1.3.12
- showtext-0.9.5
- reticulate-1.26.9


## Scripts

### Isoseq script :
- Isoseq_script/Step1_CCS.sh     : Consensus generation
- Isoseq_script/Step2_FLNC.sh    : Primer Removal
- Isoseq_script/Step3_CLUSTER.sh : Isoseq3 cluster & Polish
- Isoseq_script/Step4_CDHIT.sh   : Cluster hq_transcripts
- Isoseq_script/Step5_Annot.sh   : Gene Annotation

### Seurat script :
1. Seurat_script/Step1_DoubletFinder.sh
	- Package require :
		```
		Seurat, DoubletFinder, dplyr, ggplot2, patchwork, future
		```
	- Shell Script : 
		```shell
		Rscript script/Doublet.R data/rawdata/Sample Sample script/Seurat_lib.R data/doubletfinder/Sample
		```
	- Steps
		```
		1. Find Doublet Cells
		```

2. Seurat_script/Step2_Seurat.sh
	- Package require :
		```
		Seurat, harmony, dplyr, ggplot2, patchwork, future, Cairo, extrafont, 
		extrafontdb, ggcor, randomcoloR, Rttf2pt1, showtext
		```
	- Shell Script :
		```shell
		sed -n '/Doublet/p' data/doubletfinder/*/DF.classify.*xls | cut -f 1 > data/doubletfinder/doublet.xls
		Rscript script/Seurat.R shell/parameter.yaml data script/Seurat_lib.R
		```
	- Steps :
		```
		1. Filter Cells (remove Doublet Cells, nCount_RNA < 44000, 200 < nFeature_RNA < 9400, percent.mito < 10)
		2. Filter Data Stat
		3. Normalization Data
		4. Reduce dimension (Umap and tSNE)
		5. RunHarmony (group.by.vars = "orig.ident", project.dim = FALSE)
		6. DoFindClusters (reduction = "harmony", resolution = 0.5)
		7. Draw figures (Umap, tSNE)
		8. Find marker genes (logfc = 0.25, min_pct = 0.25, pvalue = 0.01)
		9. Plot marker genes (Heatmap, Dotplot, DensityPlot, ExpPlot, ViolinPlot)
		```

### CellChat_script
1. Step1_Seurat2CellChat.sh : 
	- Package require :
		```
		Seurat, CellChat
		```
	- Shell Script :
		```shell
		Rscript script/seurat2cellchat.R data/obj.Rda shell/cellchat.yaml CellChat_result
		```
	- Steps :
		```
		1. Read DataBase
		2. Convert Seurat object to CellChat object
		```
2. Step2_Run_CellChat.sh :
	- Package require :
		```
		Seurat, CellChat, reticulate
		```
	- Shell Script :
		```shell
		Rscript script/cellchat_sample.R CellChat_result/cellchat.init.Rds CellChat_result
		```
	- Steps :
		```
		1. identifyOverExpressedInteractions
		2. computeCommunProb
		3. aggregateNet (thresh = 0.05)
		4. computeCommunProbPathway
		5. draw figures (network, dotplot, contribution)
		```

### Seurat object Rdata file

Rda file 'obj.Rda' can be download from https://www.kdocs.cn/view/l/cdUfKEu1jIl2

