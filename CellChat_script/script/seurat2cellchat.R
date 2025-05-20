
library(Seurat)
library(CellChat)

bin <- dirname(normalizePath(sub('--file=', '',  grep('--file=', commandArgs(), value = T))))
source(paste0(bin, "/cellchat_lib.R"))
source(paste0(bin, "/modify_cellchat.R"))


args <- commandArgs(T)
if ( length(args) < 2 ) {
		stop("Rscipt foo.R <obj_file> <conf> <outdir>")
		q()
}

file <- args[1]
conf_file <- args[2]
outdir <- args[3]


dir.create(outdir, F, T)
outdir <- normalizePath(outdir)


parameter <- yaml::yaml.load_file(conf_file)

if ( is.null(parameter$sample$col_name) || is.null(parameter$cluster$col_name) ) {
		stop("\nsample$col_name and cluster$col_name must not be NULL")
		q()
}

#### find database
database <- parameter$database
if ( is.null(database) ){
		if ( ! is.null(parameter$CellChat$database) ) {
			database <- parameter$CellChat$database
		} else {
			stop("\ndatabase must not be NULL")
			q()
		}
}
if ( file.exists(database) ) {
		CellChatDB <- readRDS(database)
} else {
		CellChatDB <- switch(database,
			human = CellChatDB.human,
			mouse = CellChatDB.mouse,
			zebrafish = CellChatDB.zebrafish,
			stop("\ndefault CellChatDB only accept : human | mouse | zebrafish"))
}


#### load Seurat Object
obj <- Load(file)
print(head(obj@meta.data))


#### setting
## default assay
if ( ! is.null(parameter$data$assay) ) {
		DefaultAssay(obj) <- parameter$data$assay
}
cells.use <- Cells(obj)

## new sample name
col_name <- parameter$sample$col_name
if ( ! is.null(col_name) ) {
	if ( ! is.factor(obj@meta.data[[col_name]]) ) 
		obj@meta.data[[col_name]] <- as.factor(obj@meta.data[[col_name]])
	if ( ! is.null(parameter$sample$rename) ) {
		parameter$sample$rename <- unlist(parameter$sample$rename)
		name <- intersect(levels(obj@meta.data[[col_name]]), names(parameter$sample$rename))
		index <- levels(obj@meta.data[[col_name]]) %in% name
		levels(obj@meta.data[[col_name]])[index] <- parameter$sample$rename[name]
	}
	if ( ! is.null(parameter$sample$use)  ) {
		cells.use <- intersect(cells.use, rownames(obj@meta.data)[obj@meta.data[[parameter$sample$col_name]] %in% parameter$sample$use])
	}
}

## new cluster name
col_name <- parameter$cluster$col_name
if ( ! is.null(col_name) ) {
	if ( ! is.null(parameter$cluster$cell.set) ) {
		source('/public2/Bio/pipeline/SingleCell_Collections/SCellWare/current/R/lib/online_lib.R')
		source('/public2/Bio/pipeline/SingleCell_Collections/SCellWare/current/R/lib/utilities.R')
		obj <- CellSetFilter(obj, parameter$cluster$cell.set, col.name = col_name)
	}
	if ( ! is.factor(obj@meta.data[[col_name]]) )
		obj@meta.data[[col_name]] <- as.factor(obj@meta.data[[col_name]])
	if ( ! is.null(parameter$cluster$rename) ) {
		parameter$cluster$rename <- unlist(parameter$cluster$rename)
		name <- intersect(levels(obj@meta.data[[col_name]]), names(parameter$cluster$rename))
		index <- levels(obj@meta.data[[col_name]]) %in% name
		levels(obj@meta.data[[col_name]])[index] <- parameter$cluster$rename[name]
	}
	Idents(obj) <- col_name
	if ( ! is.null(parameter$cluster$use) ) {
		cells.use <- intersect(cells.use, rownames(obj@meta.data)[obj@meta.data[[parameter$cluster$col_name]] %in% parameter$cluster$use])
	}
}



#### filter
## cells
if ( ! is.null(parameter$cell.use) ) {
		cells <- do.call(c, lapply(parameter$cell.use, readLines))
		cells.use <- intersect(cells.use, cells)
}

## features
if ( ! is.null(parameter$features$use) ) {
		features <- readLines(parameter$features$use)
} else {
		features <- NULL
}

obj <- obj[features, cells.use]
obj@meta.data <- droplevels(obj@meta.data)



#### rename features to fit CellChatDB
features <- rownames(obj)
gene.db <- extractGene(CellChatDB)
if ( ! is.null(parameter$features$rename) ) {
		if ( file.exists(parameter$features$rename) ) {
				new_name <- read.table(parameter$features$rename, row.names = 1, stringsAsFactors = F)
		} else if ( parameter$features$rename %in% colnames(obj@misc$fdata) ) {
				new_name <- obj@misc$fdata[, parameter$features$rename, drop = F]
		} else {
				stop()
				q()
		}
		rownames(new_name) <- gsub('_', '-', rownames(new_name))

		new_name <- new_name[intersect(rownames(new_name), features), 1, drop = F]
		if ( nrow(new_name) > 0 ) {
				index <- features %in% rownames(new_name)
				features[index] <- as.character(new_name[features[index], 1])
		}		
}
names(features) <- toupper(features)
names(gene.db) <- toupper(gene.db)
index <- intersect(names(features), names(gene.db))
features[index] <- gene.db[index]

message('Number of genes in CellChatDB : ', length(gene.db))
message('Number of genes in Seurat object : ', length(features) )
message('Intersect of above : ', length(intersect(gene.db, features)))



#### split Seurat Object
#if ( is.null(parameter$split.by) ) {
#		obj.list <- list(all = obj)
#} else {
#		obj.list <- SplitObject(obj, split.by = parameter$split.by)
#}
if ( is.null(parameter$sample$group) ) {
		obj.list <- SplitObject(obj, split.by = parameter$sample$col_name)
} else if ( is.character(parameter$sample$group) ) {
		obj.list <- list(all = obj)
		names(obj.list) <- parameter$sample$group
} else if ( is.list(parameter$sample$group) ) {
		col_name <- parameter$sample$col_name
		samples <- levels(obj@meta.data[[col_name]])
		if ( '__ALLSAMPLES__' %in% names(parameter$sample$group) ) {
				parameter$sample$group[['__ALLSAMPLES__']] <- NULL
				for ( i in samples ) {
						parameter$sample$group[[i]] <- i
				}
		}
		obj.list <- list()
		for ( i in names(parameter$sample$group) ) {
				if ( is.null(parameter$sample$group[[i]]) ) {
						obj.list[[i]] <- obj
				} else {
						index <- obj@meta.data[[col_name]] %in% parameter$sample$group[[i]]
						if ( sum(index) > 0 ) {
								obj.list[[i]] <- obj[, index]
						} else {
								message("[Warning] group : '", i, "' no cells found") 
						}
				}
		}
} else {
		stop()
		q()
}


#### create CellChat Object
outfile <- list()
for ( i in names(obj.list) ) {
		message('Split object : ', i)
		outpref <- paste0(outdir, "/", i)
		dir.create(outpref, F, T)

#		cellchat <- createCellChat(object = obj.list[[i]], group.by = parameter$label.by)
		cellchat <- createCellChat(object = obj.list[[i]], group.by = parameter$cluster$col_name)
		cellchat@DB <- CellChatDB

		if ( length(features) > 0 ) {
				rownames(cellchat@data) <- features
		}
		cellchat <- subsetData(cellchat)
		cellchat@options$select_commun <- if ( 'CellChat' %in% names(parameter) ) parameter$CellChat$select_commun else parameter$select_commun

		saveRDS(cellchat, file = paste0(outpref, "/", "cellchat.init.Rds"))

		outfile[[i]] <- paste0(outpref, "/", "cellchat.init.Rds")
}
write.table(t(as.data.frame(outfile, check.names = F)), file = paste0(outdir, "/sample.list"), quote = F, col.names = F, sep = "\t")

message('ALL DONE!')


