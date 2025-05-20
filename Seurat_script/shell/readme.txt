Software :
	R-4.1.3
	DoubletFinder-2.0.3
	Seurat-4.1.0

Script :
	Step1_DoubletFinder.sh : 
		1. Find Doublet Cells

	Step2_Seurat.sh        : 
		1. Filter Cells (remove Doublet Cells, 0 < nCount_RNA < 44000, 200 < nFeature_RNA < 9400, percent.mito < 10)
		2. Filter Data Stat
		3. Normalization Data
		4. Reduce dimension (Umap and tSNE)
		5. RunHarmony (group.by.vars = "orig.ident", project.dim = FALSE)
		6, DoFindClusters (reduction = "harmony", resolution = 0.5)
		7. Draw figures (Umap, tSNE)
		8. Find marker genes (logfc = 0.25, min_pct = 0.25, pvalue = 0.01)
		9. Plot marker genes (Heatmap, Dotplot, DensityPlot, ExpPlot, ViolinPlot)

