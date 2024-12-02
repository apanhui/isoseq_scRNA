sed -n '/Doublet/p' doubletfinder/*/DF.classify.*xls | cut -f 1 > doubletfinder/doublet.xls
Rscript Seurat.R parameter.yaml Seurat Seurat_lib.R
