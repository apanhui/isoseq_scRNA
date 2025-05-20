
readRDX <- function(file) {
		con <- gzfile(file)
		on.exit(close(con))
		magic <- readChar(con, 5L, useBytes = TRUE)
		if ( grepl("RD[ABX][2-9]\n", magic) ){
				object <- get(load(file))
		} else {
				object <- readRDS(file)
		}
		return(object)
}

Load <- function(file) {
		object <- readRDX(file)
		if ( "version" %in% slotNames(object) ) {
				if ( grepl('^2', object@version) ) {
						object <- Seurat::UpdateSeuratObject(object)
				}
		}
		return(object)
}


WriteTable <- function(x, file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE, ...){
		write.table(x, file = file, quote = quote, sep = sep, row.names = row.names, col.names = col.names, ...  )
}





Scanpy2CellChat <- function(h5ad_file) {
		library(reticulate)
		ad <- import("anndata", convert = FALSE)
		ad_object <- ad$read_h5ad(h5ad_file)
		# access normalized data matrix
		data.input <- t(py_to_r(ad_object$X))
		rownames(data.input) <- rownames(py_to_r(ad_object$var))
		colnames(data.input) <- rownames(py_to_r(ad_object$obs))
		# access meta data
		meta <- py_to_r(ad_object$obs)
		cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
		return(cellchat)
}

TenX2CellChat <- function(dir = NULL, matrix_path = paste0(dir, "/filtered_feature_bc_matrix.h5"),
		cluster_path = paste0(dir, "/analysis/clustering/graphclust/clusters.csv") ) {
		source("/public2/Bio/pipeline/SingleCell_Collections/SCellWare/current/R/lib/data_import.R")
		data <- .Read10X(matrix_path, use.names = TRUE)
		data <- normalizeData(data)
		sep <- if ( grepl('.csv$', cluster_path) ) ',' else '\t'
		meta <- read.table(cluster_path, header = T, row.names = 1, sep = sep)		
		if ( all(grepl(pattern = "-[0-9]+$", x = rownames(meta))) ){
				rownames(meta) <- as.vector(x = as.character(x = sapply(X = rownames(meta), FUN = Seurat:::ExtractField, field = 1, delim = "-")))
		}
		cellchat <- createCellChat(object = data, meta = meta, group.by = colnames(meta)[1])
		return(cellchat)
}

StatCommunProb <- function(object, outpref = "LR.CommunProb", slot.name = 'net', color = NULL, ...) {
		data <- subsetCommunication(object, min.prob = -1, slot.name = slot.name, ...)
		subdata <- subsetCommunication(object, min.prob = 0, slot.name = slot.name, ...)
		if ( object@options$mode == 'single' ) {
				data$source <- factor(data$source, levels = levels(object@idents))
				data$target <- factor(data$target, levels = levels(object@idents))
				subdata$source <- factor(subdata$source, levels = levels(object@idents))
				subdata$target <- factor(subdata$target, levels = levels(object@idents))

				.StatCommunProb(data, subdata, outpref = outpref, slot.name = slot.name, color = color)

				pdf(paste0(outpref, ".Heatmap.pdf"), 12, 6)
				ht1 <- netVisual_heatmap(object, color.heatmap = c("white", "darkred"))
				ht2 <- netVisual_heatmap(object, color.heatmap = c("white", "darkblue"), measure = "weight")
				ComplexHeatmap::draw(ht1 + ht2, ht_gap = unit(0.5, "cm"), legend_gap = grid::unit(3, 'line'))
			dev.off()
		} else if ( object@options$mode == 'merged' ) {
				for ( i in names(data) ) {
						data[[i]]$source <- factor(data[[i]]$source, levels = levels(object@idents))
						data[[i]]$target <- factor(data[[i]]$target, levels = levels(object@idents))
						subdata[[i]]$source <- factor(subdata[[i]]$source, levels = levels(object@idents))
						subdata[[i]]$target <- factor(subdata[[i]]$target, levels = levels(object@idents))

						.StatCommunProb(data[[i]], subdata[[i]], sample = i, outpref = outpref, slot.name = slot.name, color = color)
				}
		} else {
				stop()
		}
		invisible(subdata)
}

.StatCommunProb <- function(data, subdata, sample = NULL, outpref = "LR.CommunProb", slot.name = 'net', color = NULL, top = 15) {
	if ( ! is.null(sample) ) outpref <- paste0(outpref, ".", sample)
	if ( is.null(color) ) color <- c("lightgrey", "red")

	if ( slot.name == 'net' ) {
			data <- data %>% select(-interaction_name, -evidence) %>% rename(Pair = interaction_name_2)
			subdata <- subdata %>% select(-interaction_name, -evidence) %>% rename(Pair = interaction_name_2)
			stat <- subset(data, prob > 0) %>% group_by(source, target) %>% summarise(count = n(), prob = sum(prob))
			write.table(stat, file = paste0(outpref, ".stat.xls"), row.names = F, sep = "\t", quote = F)
			p <- ggplot(stat, aes(x = source, y = target)) +
				geom_point(aes(color = prob, size = count)) +
				scale_color_gradientn(colours = color) +
				theme_classic() +
				theme(panel.border = element_rect(fill = NA), axis.line = element_blank(),
						axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
			ggsave(p, file = paste0(outpref, ".stat.dotplot.pdf"), width = 6, height = 6)

	} else if ( slot.name == 'netP' ) {
#			subdata <- subset(data, prob > 0)
			stat <- rbind(
					subdata %>% group_by(Cluster = source, pathway_name) %>% summarise(prob = sum(prob)) %>% mutate(type = 'source'),
					subdata %>% group_by(Cluster = target, pathway_name) %>% summarise(prob = sum(prob)) %>% mutate(type = 'target')
				) %>% select(Cluster, type, pathway_name, prob)
			write.table(stat, file = paste0(outpref, ".stat.xls"), row.names = F, sep = "\t", quote = F)			

			df <- subdata %>% group_by(pathway_name) %>%
				summarise(Pval = sum(pval), Prob = sum(prob)) %>% arrange(Pval, -Prob) %>%
				head(top) %>% left_join(y = subdata) %>%
				group_by(pathway_name) %>% top_n(top, prob)
			df$name <- paste(df$source, df$target, sep = ' -> ')

			p <- ggplot(df, aes(x = pathway_name, y = name)) +
			    geom_point(aes(color = prob, size = -pval)) +
				scale_color_gradientn(colours = color) +
				labs(x = 'Pathway', y = 'source-target') + 
				theme_classic() +
				theme(panel.border = element_rect(fill = NA), axis.line = element_blank(),
						axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
			p <- p + scale_size("pval",
						labels = unique(rev(scales::breaks_extended()(range(df$pval)))),
						breaks = unique(scales::breaks_extended()(range(-df$pval))))			
			p <- p + coord_flip()
			height <- max(6, length(unique(df$pathway_name)) * 0.2)
			width <- max(6, length(unique(df$name)) * 0.2)
			ggsave(p, file = paste0(outpref, ".stat.dotplot.pdf"), width = width, height = height, limitsize = FALSE)
	} else {
			stop()
	}
	write.table(data, file = paste0(outpref, ".all.xls"), row.names = F, sep = "\t", quote = F)
	write.table(subdata, file = paste0(outpref, ".xls"), row.names = F, sep = "\t", quote = F)


}


.PlotNet <- function(data, sample = NULL, split.group = FALSE, outpref = "LR.interaction", groupSize = 20, is.combine = FALSE, ...) {
		if ( ! is.null(sample) ) outpref <- paste0(outpref, ".", sample)

		.PlotNet2(data$count, type = 'count', outpref = outpref, split.group = split.group, is.combine = is.combine, groupSize = groupSize, ...)

		.PlotNet2(data$weight, type = 'weight', outpref = outpref, split.group = split.group, is.combine = is.combine, groupSize = groupSize, ...)
}


.PlotNet2 <- function(data, type = 'count', split.group = FALSE, outpref = "LR.interaction", groupSize = 20, is.combine = FALSE, ...) {
		title.name <- if ( type == 'count' ) "Numberof interactions" else " Interaction weights/strength"

		point.num <- length(unique(rownames(data), colnames(data)))
		w <- h <- max(6, point.num * 3 / 20 + 3)
		if ( ! split.group ) {
				pdf(paste0(outpref, ".", type, ".pdf"), width = w, height = h)
				netVisual_circle(data, vertex.weight = groupSize, weight.scale = T, title.name = title.name, ...)
				dev.off()
		} else {
				plots <- list()
				if ( is.combine ) {
						ncol <- ceiling(sqrt(nrow(data)))
						nrow <- ceiling(nrow(data)/ncol)
						pdf(paste0(outpref, ".InCluster.", type, ".pdf"), width = 4/6 * w * ncol, height = 4/6 * h * nrow)
						par(mfrow = c(nrow, ncol), xpd = TRUE)
				}
				for (i in 1:nrow(data)) {
						mat <- matrix(0, nrow = nrow(data), ncol = ncol(data), dimnames = dimnames(data))
						mat[i, ] <- data[i, ]
						if ( ! is.combine ) {
								pdf(paste0(outpref, ".", rownames(data)[i], ".", type, ".pdf"), width = w, height = h)
								title.name <- rownames(data)[i]
						}
						netVisual_circle(mat, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(data), title.name = title.name, ...)
						if ( ! is.combine ) {
								dev.off()
						}
				}
				if ( is.combine ) {
						dev.off()
				}
		}
}


PlotNet <-function(object, outpref = "LR.interaction", ...) {
		if ( object@options$mode == 'single' ) {
				groupSize <- as.numeric(table(object@idents))
				.PlotNet(object@net, outpref = outpref, groupSize = groupSize, ...)
		} else if ( object@options$mode == 'merged' ) {
				for ( i in names(object@net) ) {
						groupSize <- as.numeric(table(object@idents[[i]]))
						.PlotNet(object@net[[i]], sample = i, outpref = outpref, groupSize = groupSize, ...)
				}
		} else {
				stop()
		}
}


PlotLRDotPlot <- function(object, outpref = "LR.DotPlot", top = 5, thresh = 0.05, ... ) {
		data <- subsetCommunication(object, thresh = thresh, ...)
		if ( !is.null(top) && top > 0 ) {
			sources.use <- (data %>% group_by(source) %>% summarise(PROB = sum(prob)) %>% arrange(-PROB) %>% head(top))[[1]]
			targets.use <- (data %>% group_by(target) %>% summarise(PROB = sum(prob)) %>% arrange(-PROB) %>% head(top))[[1]]
			pairLR.use <- (subset(data, source %in% sources.use & target %in% targets.use) %>% group_by(source, target) %>% arrange(-prob) %>% filter(1:n() <= top))[, 'interaction_name', drop = FALSE]
			outpref <- paste0(outpref, ".top", top)
		} else {
			sources.use <- targets.use <- pairLR.use <- NULL
		}
		p <- netVisual_bubble(object, pairLR.use = pairLR.use, sources.use = sources.use, targets.use = targets.use, remove.isolate = T, thresh = thresh)
		p <- p + coord_equal()
#		p$data <- p$data[p$data$source.target %in% paste(stat$source, stat$target, sep = ' -> '), ]
#		p$data <- droplevels(p$data)
		color.use <- rev(RColorBrewer::brewer.pal(n = 10, name = 'Spectral'))
		if (min(p$data$prob, na.rm = T) != max(p$data$prob, na.rm = T)) {
				p <- p + scale_colour_gradientn(colors = colorRampPalette(color.use)(99), 
						limits = c(quantile(p$data$prob, 0, na.rm = T), quantile(p$data$prob, 1, na.rm = T)), 
						breaks = c(quantile(p$data$prob, 0, na.rm = T), quantile(p$data$prob, 1, na.rm = T)),
						na.value = "white", labels = c("min", "max")) +
					guides(color = guide_colourbar(barwidth = 0.5, title = "Commun. Prob."))
		}
		else {
				p <- p + scale_colour_gradientn(colors = colorRampPalette(color.use)(99), na.value = "white") +
					guides(color = guide_colourbar(barwidth = 0.5, title = "Commun. Prob."))
		}
		width <- max(7, nlevels(p$data$source.target) * 0.2 + 1)
		height <- max(6, nlevels(p$data$interaction_name_2) * 0.2)
		ggsave(p, file = paste0(outpref, "_prob.pdf"), width = width, height = height, limitsize = FALSE)
}

GetWH <- function(gg, scale = 0.12, scale.width = scale, scale.height = scale,
			default.size = 6, default.width = default.size, default.height = default.size) {
		x <- as.character(gg$mapping$x)[2]
		y <- as.character(gg$mapping$y)[2]
		width <- max(default.width, length(unique(gg$data[[x]])) * scale.width)
		height <- max(default.height, length(unique(gg$data[[y]])) * scale.height)
		return(c(width, height))
}

PlotDiffCommunProb <- function(cellchat, outpref = "Diff.CommunProb", color.text = c("#377EB8", "#FF7F00"), ...) {
	dt <- netVisual_bubble(cellchat, comparison = c(1,2), angle.x = 45, remove.isolate = T, color.text = color.text, return.data = TRUE, ...)
	p <- dt[['gg.obj']]
	wh <- GetWH(p)
	ggsave(p, file = paste0(outpref, ".bubble.pdf"), width = wh[1] + 1.5, height = wh[2] + 1.5, limitsize = FALSE)

	for ( i in c(2) ) {
		name <- names(cellchat@net)[i]

		j <- switch(i, '1' = 2, '2' = 1)
		gg1 <- netVisual_bubble(cellchat, comparison = c(1,2), max.dataset = i, title.name = paste0("Increased signaling in ", name), angle.x = 45, remove.isolate = T, color.text = color.text, ...)
		wh <- GetWH(gg1)
		ggsave(gg1, file = paste0(outpref, ".Increased.bubble.pdf"), width = wh[1] + 1.5, height = wh[2] + 1.5, limitsize = FALSE)

		gg2 <- netVisual_bubble(cellchat, comparison = c(1,2), max.dataset = j, title.name = paste0("Decreased signaling in ", name), angle.x = 45, remove.isolate = T, color.text = color.text, ...)
		wh <- GetWH(gg2)
		ggsave(gg2, file = paste0(outpref, ".Decreased.bubble.pdf"), width = wh[1] + 1.5, height = wh[2] + 1.5, limitsize = FALSE)
	}

	invisible(dt[['communication']])
}

PlotDiffExp <- function(cellchat, pairLR.use = NULL, outpref = "Diff.Exp", ...) {
	p <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, comparison = c(1, 2),  angle.x = 45, remove.isolate = T,title.name = "Up-regulated signaling", ...)
	wh <- GetWH(p)
	ggsave(p, file = paste0(outpref, ".bubble.pdf"), width = wh[1], height = wh[2], limitsize = F)
}


StatDiffNumber_xls <- function(cellchat,
			sample1 = names(cellchat@net)[1], sample2 = names(cellchat@net)[2],
			net1 = cellchat@net[[sample1]], net2 = cellchat@net[[sample2]],
			count1 = net1$count, count2 = net2$count,
			weight1 = net1$weight, weight2 = net2$weight,
			outpref = NULL) {
		dt <- reshape2::melt(count1, varnames = c('Source', 'Target'), value.name = paste0('Count_', sample1)) %>%
			full_join(y = reshape2::melt(count2, varnames = c('Source', 'Target'), value.name = paste0('Count_', sample2)))
		dt[[paste0('Count_', sample1, '-vs-', sample2)]] <- dt[[4]] - dt[[3]]
		dt <- dt %>% full_join(y = reshape2::melt(weight1, varnames = c('Source', 'Target'), value.name = paste0('Prob_', sample1))) %>%
			full_join(y = reshape2::melt(weight2, varnames = c('Source', 'Target'), value.name = paste0('Prob_', sample2)))
		dt[[paste0('Prob_', sample1, '-vs-', sample2)]] <- dt[[7]] - dt[[6]]
		if ( ! is.null(outpref) ) {
				WriteTable(dt, paste0(outpref, ".xls"))
		}
		invisible(dt)
}

StatDiffNumber <- function(cellchat, outpref = 'Diff.Number') {
		dt <- StatDiffNumber_xls(cellchat, outpref = 'Diff.Number')

		dt <- dt[, c(1,2,5,8)]
		colnames(dt) <- c('Source', 'Target', 'diff_count', 'diff_prob')
		#p <- ggplot(dt, aes(x = Source, y = Target, size = diff_count, color = diff_prob)) +
		#			geom_point()
		# TODO

		gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
		gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
		ggsave(gg1 + gg2, file = paste0(outpref, ".barplot.pdf"), width = 12, height = 6)

		pdf(paste0(outpref, ".network.pdf"), 12, 6)
		par(mfrow = c(1,2), xpd=TRUE)
		netVisual_diffInteraction(cellchat, weight.scale = T)
		netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
		dev.off()

		pdf(paste0(outpref, ".heatmap.pdf"), 12, 6)
		ht1 <- netVisual_heatmap(cellchat, legend.name = 'Relative Number')
		ht2 <- netVisual_heatmap(cellchat, measure = "weight", legend.name = 'Relative strength')
		ComplexHeatmap::draw(ht1 + ht2, ht_gap = unit(0.5, "cm"), legend_gap = grid::unit(3, 'line'))
		dev.off()

		invisible(dt)
}

GetCentr <- function(cellchat, sample = NULL, netP = NULL) {
	if ( is.null(netP) ) {
			if ( cellchat@options$mode == 'merged' ) {
					if ( is.null(sample) ) {
						stop("when 'netP' is NULL and cellchat@options$mode == 'merged', 'sample' must not be NULL")
					}
					netP <- cellchat@netP[[sample]]
			} else {
					netP <- cellchat@netP
			}
	}
	if ( is.null(netP$centr) ) {
			stop("you should run netAnalysis_computeCentrality first")
	}
	dt <- do.call(rbind,
		lapply(names(netP$centr),
			function(j) {
				x <- netP$centr[[j]]
				y <- data.frame(pathway = j, outgoing = x$outdeg, incoming = x$indeg)
				y$cluster <- rownames(y)
				return(y)
			}
		)
	)
	return(dt)
}

PlotPathwayDiff_network <- function(cellchat, pathways = NULL, outpref = "Diff.Pathway", plot.size = 6) {
	if ( is.null(pathways) ) {
			pathways <- union(cellchat@netP[[1]]$pathways, cellchat@netP[[2]]$pathways)
	}
	for ( i in pathways ) {
		if ( i %in% cellchat@netP[[1]]$pathways ) {
				net1 <- cellchat@netP[[1]]$prob[,,i]
		} else {
				net1 <- matrix(0, nrow = nlevels(cellchat@idents$joint), ncol = nlevels(cellchat@idents$joint))
		}
		if ( i %in% cellchat@netP[[2]]$pathways ) {
				net2 <- cellchat@netP[[2]]$prob[,,i]
		} else {
				net2 <- matrix(0, nrow = nlevels(cellchat@idents$joint), ncol = nlevels(cellchat@idents$joint))
		}
		net.diff <- net1 - net2

		pdf(paste(outpref, i, "network.pdf", sep = "."), plot.size, plot.size)
		netVisual_diffInteraction(net.diff = net.diff, weight.scale = T, measure = "weight")
		dev.off()
	}
}

PlotSignalingPatterns_heatmap <- function(object.list, pathways = NULL, outpref = 'SignalingPatterns', idents = NULL) {
	if ( is.null(idents) ) {
		idents <- union(levels(object.list[[1]]@idents), levels(object.list[[2]]@idents))
	} else {
		object.list[[1]] <- subsetCellChat(object.list[[1]], idents = idents)
		object.list[[2]] <- subsetCellChat(object.list[[2]], idents = idents)
	}
	if ( is.null(pathways) ) {
		pathways <- union(object.list[[1]]@netP$pathways, object.list[[2]]@netP$pathways)
	}
	width  <- max(7, length(idents) * 0.12) ## in inch
	height <- max(7, length(pathways) * 0.12) 
	for ( i in c("outgoing", "incoming", "all") ) {
		color.heatmap <- switch(i, outgoing = "BuGn", incoming = "GnBu", all = "OrRd")
		ht1 <- netAnalysis_signalingRole_heatmap(object.list[[1]], pattern = i, signaling = pathways, title = names(object.list)[1], width = width * 2.54, height = height * 2.54, color.heatmap = color.heatmap) ## in cm
		ht2 <- netAnalysis_signalingRole_heatmap(object.list[[2]], pattern = i, signaling = pathways, title = names(object.list)[2], width = width * 2.54, height = height * 2.54, color.heatmap = color.heatmap)
		pdf(paste0(outpref, '.', i ,".pdf"), (width + 2 ) * 2, height + 2)
		ComplexHeatmap::draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
		dev.off()
	}
}

StatSignalingPatterns <- function(cellchat, outpref = 'SignalingPatterns') {
	dt <- list()
	for ( i in samples) {
		dt[[i]] <- GetCentr(cellchat, i)[, c(1,4,2,3)]
		dt[[i]]$all <- dt[[i]][[3]] + dt[[i]][[4]]
		name <- paste0("Prob_", i)
		colnames(dt[[i]]) <- c("Pathway", "Cluster", paste0(c('Outgoing', 'Incoming', 'Overall'), name))
}
	dt <- full_join(dt[[1]], dt[[2]])
	dt[is.na(dt)] <- 0
	WriteTable(dt, file = paste0(outpref, ".xls"))
}


SelectCommunication <- function(cellchat, source = NULL, target = NULL,
			method = c('both_side', 'one_side')) {
	method <- match.arg(method)
	if ( method == 'both_side' ){
		if ( !is.null(source) ) {
			source <- intersect(source, levels(cellchat@idents))
			cellchat@net$prob[source, source, ] <- 0
		}
		if ( !is.null(target) ) {
			target <- intersect(target, levels(cellchat@idents))
			cellchat@net$prob[target, target, ] <- 0
		}
	} else if ( method == 'one_side' ) {
		if ( is.null(source) || is.null(target) ) {
			return(cellchat)
		}
		source <- intersect(source, levels(cellchat@idents))
		target <- intersect(target, levels(cellchat@idents))
		non_source <- setdiff(levels(cellchat@idents), source)
		non_target <- setdiff(levels(cellchat@idents), target)
		cellchat@net$prob[non_source, , ] <- 0
		cellchat@net$prob[, non_target, ] <- 0
	}
	return(cellchat)
}



CalLRAvg <- function(object) {
	data <- as.matrix(object@data.signaling)
	pairLRsig <- object@LR$LRsig
	complex_input <- object@DB$complex
	cofactor_input <- object@DB$cofactor
	group <- object@idents

	geneL <- as.character(pairLRsig$ligand)
	geneR <- as.character(pairLRsig$receptor)

	data.use <- data/max(data)
	data.use.avg <- aggregate(t(data.use), list(group), FUN = CellChat::triMean)
	data.use.avg <- t(data.use.avg[, -1])
	colnames(data.use.avg) <- levels(group)

	dataLavg <- computeExpr_LR(geneL, data.use.avg, complex_input)
	dataRavg <- computeExpr_LR(geneR, data.use.avg, complex_input)
	dataRavg.co.A.receptor <- computeExpr_coreceptor(cofactor_input,
		data.use.avg, pairLRsig, type = "A")
	dataRavg.co.I.receptor <- computeExpr_coreceptor(cofactor_input,
		data.use.avg, pairLRsig, type = "I")
	dataRavg <- dataRavg * dataRavg.co.A.receptor/dataRavg.co.I.receptor

	rownames(dataLavg) <- geneL
	rownames(dataRavg) <- geneR
	colnames(dataLavg) <- levels(group)
	colnames(dataRavg) <- levels(group)
	return(list(
		ligand = dataLavg[unique(geneL),,drop = F],
		receptor = dataRavg[unique(geneR),,drop = F]
		))
}

GetOnlineData <- function(cellchat) {
	data <- subsetCommunication(cellchat, min.prob = -1, slot.name = 'net')
	data <- data %>% select(-interaction_name, -evidence) %>%
		select(-pathway_name, -annotation) %>%
		rename(Pair = interaction_name_2, Probability = prob, Pvalue = pval) %>%
		rename(Source = source, Target = target, Ligand = ligand, Receptor = receptor)

	avg <- CalLRAvg(cellchat)
	avgL <- reshape2::melt(avg$ligand, varnames = c('Ligand', 'Source'), value.name = 'Source_mean_exp')
	avgR <- reshape2::melt(avg$receptor, varnames = c('Receptor', 'Target'), value.name = 'Target_mean_exp')
	avgL$Ligand   <- as.character(avgL$Ligand)
	avgL$Source   <- as.character(avgL$Source)
	avgR$Receptor <- as.character(avgR$Receptor)
	avgR$Target   <- as.character(avgR$Target)

	data <- data %>% left_join(y = avgL) %>% left_join(y = avgR) %>%
		arrange(desc(Source_mean_exp * Target_mean_exp)) %>%
		select(Source, Target, Ligand, Receptor, Pair, Source_mean_exp, Target_mean_exp, Probability, Pvalue)
	return(data)
}


SubsetObj <- function(object, cells = NULL, sample = NULL, cluster = NULL, sample.name = "orig.ident", cluster.name = "seurat_clusters", ...) {
		cells <- if ( ! is.null(cells) ) cells else Cells(object)
		if ( ! is.null(sample) ) {
				if ( is.list(sample) ) {
				## TODO : rename sample
				} else {
						cells.sample <- Cells(object)[object[[sample.name]][[1]] %in% sample]
				}
				cells <- intersect(cells, cells.sample)
		}
		if ( ! is.null(cluster) ) {
				if ( is.list(cluster) ) {
				## TODO : rename cluster
				} else {
						cells.cluster <- Cells(object)[object[[cluster.name]][[1]] %in% cluster]
				}
				cells <- intersect(cells, cells.cluster)
		}
		others <- list(...)
		if ( length(others) > 0 ) {

		}
		object <- object[, cells]
		object@meta.data <- droplevels(object@meta.data)
		if ( 'images' %in% slotNames(object) ) {
				## when object[, cells], if image's name is like 'A-B', it will be duplicated with name 'A.B',
				## here substracting right name images with 'sample.name' in object@meta.data
				object@images <- object@images[Images(object) %in% unique(object@meta.data[[sample.name]])]
				for ( image in Images(object) ) {
						## object[, cells] did not apply to object@images
						image.cells <- intersect(rownames(object@images[[image]]@coordinates), cells)
						object@images[[image]]@coordinates <- object@images[[image]]@coordinates[image.cells, ]
				}
		}
		return(object)
}

