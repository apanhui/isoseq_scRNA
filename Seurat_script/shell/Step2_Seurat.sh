sed -n '/Doublet/p' ../data/doubletfinder/*/DF.classify.*xls | cut -f 1 > ../data/doubletfinder/doublet.xls
Rscript ../script/Seurat.R parameter.yaml ../data ../script/Seurat_lib.R
