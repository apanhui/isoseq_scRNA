
setIdent <- function (object, ident.use = NULL, levels = NULL, display.warning = TRUE) 
{
    if (!is.null(ident.use)) {
        object@idents <- as.factor(object@meta[[ident.use]])
    }
    if (!is.null(levels)) {
        object@idents <- factor(object@idents, levels = levels)
    }
## here
#    if ("0" %in% as.character(object@idents)) {
#        stop("Cell labels cannot contain `0`! ")
#    }
    if (length(object@net) > 0) {
        if (all(dimnames(object@net$prob)[[1]] %in% levels(object@idents))) {
            message("Reorder cell groups! ")
            cat("The cell group order before reordering is ", 
                dimnames(object@net$prob)[[1]], "\n")
            idx <- match(levels(object@idents), dimnames(object@net$prob)[[1]])
            object@net$prob <- object@net$prob[idx, , ]
            object@net$prob <- object@net$prob[, idx, ]
            object@net$pval <- object@net$pval[idx, , ]
            object@net$pval <- object@net$pval[, idx, ]
            cat("The cell group order after reordering is ", 
                dimnames(object@net$prob)[[1]], "\n")
        }
        else {
            message("Rename cell groups but do not change the order! ")
            cat("The cell group order before renaming is ", dimnames(object@net$prob)[[1]], 
                "\n")
            dimnames(object@net$prob) <- list(levels(object@idents), 
                levels(object@idents), dimnames(object@net$prob)[[3]])
            dimnames(object@net$pval) <- dimnames(object@net$prob)
            cat("The cell group order after renaming is ", dimnames(object@net$prob)[[1]], 
                "\n")
        }
        if (display.warning) {
            warning("All the calculations after `computeCommunProb` should be re-run!!\n    These include but not limited to `computeCommunProbPathway`,`aggregateNet`, and `netAnalysis_computeCentrality`.")
        }
    }
    return(object)
}
environment(setIdent) <- asNamespace('CellChat')
assignInNamespace("setIdent", setIdent, ns = "CellChat")
#rm(setIdent)

subsetCommunication_internal <- function (net, LR, cells.level, slot.name = "net", sources.use = NULL, 
    targets.use = NULL, signaling = NULL, pairLR.use = NULL, 
    thresh = 0.05, datasets = NULL, ligand.pvalues = NULL, ligand.logFC = NULL, 
    ligand.pct.1 = NULL, ligand.pct.2 = NULL, receptor.pvalues = NULL, 
    receptor.logFC = NULL, receptor.pct.1 = NULL, receptor.pct.2 = NULL,
	min.prob = 0 ## here
	) 
{
    if (!is.data.frame(net)) {
        prob <- net$prob
        pval <- net$pval
        prob[pval > thresh] <- 0 # in <computeCommunProbPathway> : prob[pval > thresh] <- 0
        net <- reshape2::melt(prob, value.name = "prob")
        colnames(net)[1:3] <- c("source", "target", "interaction_name")
		net$source <- as.character(net$source) #
		net$target <- as.character(net$target) # here
        net.pval <- reshape2::melt(pval, value.name = "pval")
        net$pval <- net.pval$pval
        net <- subset(net, prob > min.prob)
    }
    if (!("ligand" %in% colnames(net))) {
        pairLR <- dplyr::select(LR, c("interaction_name_2", "pathway_name", 
            "ligand", "receptor", "annotation", "evidence"))
        idx <- match(net$interaction_name, rownames(pairLR))
        net <- cbind(net, pairLR[idx, ])
    }
    if (!is.null(signaling)) {
        pairLR.use <- data.frame()
        for (i in 1:length(signaling)) {
            pairLR.use.i <- searchPair(signaling = signaling[i], 
                pairLR.use = LR, key = "pathway_name", matching.exact = T, 
                pair.only = T)
            pairLR.use <- rbind(pairLR.use, pairLR.use.i)
        }
    }
    if (!is.null(pairLR.use)) {
        net <- tryCatch({
            subset(net, interaction_name %in% pairLR.use$interaction_name)
        }, error = function(e) {
            subset(net, pathway_name %in% pairLR.use$pathway_name)
        })
    }
    if (!is.null(datasets)) {
        if (!("datasets" %in% colnames(net))) {
            stop("Please run `identifyOverExpressedGenes` and `netMappingDEG` before selecting 'datasets'")
        }
        net <- net[net$datasets %in% datasets, , drop = FALSE]
    }
    if (!is.null(ligand.pvalues)) {
        if (!("ligand.pvalues" %in% colnames(net))) {
            stop("Please run `identifyOverExpressedGenes` and `netMappingDEG` before using the threshold 'ligand.pvalues'")
        }
        net <- net[net$ligand.pvalues <= ligand.pvalues, , drop = FALSE]
    }
    if (!is.null(ligand.logFC)) {
        if (!("ligand.logFC" %in% colnames(net))) {
            stop("Please run `identifyOverExpressedGenes` and `netMappingDEG` before using the threshold 'ligand.logFC'")
        }
        if (ligand.logFC >= 0) {
            net <- net[net$ligand.logFC >= ligand.logFC, , drop = FALSE]
        }
        else {
            net <- net[net$ligand.logFC <= ligand.logFC, , drop = FALSE]
        }
    }
    if (!is.null(ligand.pct.1)) {
        if (!("ligand.pct.1" %in% colnames(net))) {
            stop("Please run `identifyOverExpressedGenes` and `netMappingDEG` before using the threshold 'ligand.pct.1'")
        }
        net <- net[net$ligand.pct.1 >= ligand.pct.1, , drop = FALSE]
    }
    if (!is.null(ligand.pct.2)) {
        if (!("ligand.pct.2" %in% colnames(net))) {
            stop("Please run `identifyOverExpressedGenes` and `netMappingDEG` before using the threshold 'ligand.pct.2'")
        }
        net <- net[net$ligand.pct.2 >= ligand.pct.2, , drop = FALSE]
    }
    if (!is.null(receptor.pvalues)) {
        if (!("receptor.pvalues" %in% colnames(net))) {
            stop("Please run `identifyOverExpressedGenes` and `netMappingDEG` before using the threshold 'receptor.pvalues'")
        }
        net <- net[net$receptor.pvalues <= receptor.pvalues, 
            , drop = FALSE]
    }
    if (!is.null(receptor.logFC)) {
        if (!("receptor.logFC" %in% colnames(net))) {
            stop("Please run `identifyOverExpressedGenes` and `netMappingDEG` before using the threshold 'receptor.logFC'")
        }
        if (receptor.logFC >= 0) {
            net <- net[net$receptor.logFC >= receptor.logFC, 
                , drop = FALSE]
        }
        else {
            net <- net[net$receptor.logFC <= receptor.logFC, 
                , drop = FALSE]
        }
    }
    if (!is.null(receptor.pct.1)) {
        if (!("receptor.pct.1" %in% colnames(net))) {
            stop("Please run `identifyOverExpressedGenes` and `netMappingDEG` before using the threshold 'receptor.pct.1'")
        }
        net <- net[net$receptor.pct.1 >= receptor.pct.1, , drop = FALSE]
    }
    if (!is.null(receptor.pct.2)) {
        if (!("receptor.pct.2" %in% colnames(net))) {
            stop("Please run `identifyOverExpressedGenes` and `netMappingDEG` before using the threshold 'receptor.pct.2'")
        }
        net <- net[net$receptor.pct.2 >= receptor.pct.2, , drop = FALSE]
    }
    net <- net[rowSums(is.na(net)) != ncol(net), , drop = FALSE]
    if (nrow(net) == 0) {
        stop("No significant signaling interactions are inferred based on the input!")
    }
    if (slot.name == "netP") {
        net <- dplyr::select(net, c("source", "target", "pathway_name", 
            "prob", "pval", "annotation"))
        net$source_target <- paste(net$source, net$target, sep = "sourceTotarget")
        net.pval <- net %>% group_by(source_target, pathway_name) %>% 
            summarize(pval = mean(pval), .groups = "drop")
        net <- net %>% group_by(source_target, pathway_name) %>% 
            summarize(prob = sum(prob), .groups = "drop")
        a <- stringr::str_split(net$source_target, "sourceTotarget", 
            simplify = T)
        net$source <- as.character(a[, 1])
        net$target <- as.character(a[, 2])
        net <- dplyr::select(net, -source_target)
        net$pval <- net.pval$pval
    }
    if (!is.null(sources.use)) {
        if (is.numeric(sources.use)) {
            sources.use <- cells.level[sources.use]
        }
        net <- subset(net, source %in% sources.use)
    }
    if (!is.null(targets.use)) {
        if (is.numeric(targets.use)) {
            targets.use <- cells.level[targets.use]
        }
        net <- subset(net, target %in% targets.use)
    }
    net <- BiocGenerics::as.data.frame(net, stringsAsFactors = FALSE)
    if (nrow(net) == 0) {
        warning("No significant signaling interactions are inferred!")
    }
    else {
        rownames(net) <- 1:nrow(net)
    }
    if (slot.name == "net") {
        if (("ligand.logFC" %in% colnames(net)) & ("datasets" %in% 
            colnames(net))) {
            net <- net[, c("source", "target", "ligand", "receptor", 
                "prob", "pval", "interaction_name", "interaction_name_2", 
                "pathway_name", "annotation", "evidence", "datasets", 
                "ligand.logFC", "ligand.pct.1", "ligand.pct.2", 
                "ligand.pvalues", "receptor.logFC", "receptor.pct.1", 
                "receptor.pct.2", "receptor.pvalues")]
        }
        else if ("ligand.logFC" %in% colnames(net)) {
            net <- net[, c("source", "target", "ligand", "receptor", 
                "prob", "pval", "interaction_name", "interaction_name_2", 
                "pathway_name", "annotation", "evidence", "ligand.logFC", 
                "ligand.pct.1", "ligand.pct.2", "ligand.pvalues", 
                "receptor.logFC", "receptor.pct.1", "receptor.pct.2", 
                "receptor.pvalues")]
        }
        else {
            net <- net[, c("source", "target", "ligand", "receptor", 
                "prob", "pval", "interaction_name", "interaction_name_2", 
                "pathway_name", "annotation", "evidence")]
        }
    }
    else if (slot.name == "netP") {
        net <- net[, c("source", "target", "pathway_name", "prob", 
            "pval")]
    }
    return(net)
}
environment(subsetCommunication_internal) <- asNamespace('CellChat')
assignInNamespace("subsetCommunication_internal", subsetCommunication_internal, ns = "CellChat")
#rm(subsetCommunication_internal)

subsetCommunication <- function (object = NULL, net = NULL, slot.name = "net", sources.use = NULL, 
    targets.use = NULL, signaling = NULL, pairLR.use = NULL, 
    thresh = 0.05, datasets = NULL, ligand.pvalues = NULL, ligand.logFC = NULL, 
    ligand.pct.1 = NULL, ligand.pct.2 = NULL, receptor.pvalues = NULL, 
    receptor.logFC = NULL, receptor.pct.1 = NULL, receptor.pct.2 = NULL, ...) 
{
    if (!is.null(pairLR.use)) {
        if (!is.data.frame(pairLR.use)) {
            stop("pairLR.use should be a data frame with a signle column named either 'interaction_name' or 'pathway_name' ")
        }
        else if ("pathway_name" %in% colnames(pairLR.use)) {
            message("slot.name is set to be 'netP' when pairLR.use contains signaling pathways")
            slot.name = "netP"
        }
    }
    if (!is.null(pairLR.use) & !is.null(signaling)) {
        stop("Please do not assign values to 'signaling' when using 'pairLR.use'")
    }
    if (object@options$mode == "single") {
        if (is.null(net)) {
            net <- slot(object, "net")
        }
        LR <- object@LR$LRsig
        cells.level <- levels(object@idents)
        df.net <- subsetCommunication_internal(net, LR, cells.level, 
            slot.name = slot.name, sources.use = sources.use, 
            targets.use = targets.use, signaling = signaling, 
            pairLR.use = pairLR.use, thresh = thresh, datasets = datasets, 
            ligand.pvalues = ligand.pvalues, ligand.logFC = ligand.logFC, 
            ligand.pct.1 = ligand.pct.1, ligand.pct.2 = ligand.pct.2, 
            receptor.pvalues = receptor.pvalues, receptor.logFC = receptor.logFC, 
            receptor.pct.1 = receptor.pct.1, receptor.pct.2 = receptor.pct.2, ...)
    }
    else if (object@options$mode == "merged") {
        if (is.null(net)) {
            net0 <- slot(object, "net")
            df.net <- vector("list", length(net0))
            names(df.net) <- names(net0)
            for (i in 1:length(net0)) {
                net <- net0[[i]]
                LR <- object@LR[[i]]$LRsig
                cells.level <- levels(object@idents[[i]])
                df.net[[i]] <- subsetCommunication_internal(net, 
                  LR, cells.level, slot.name = slot.name, sources.use = sources.use, 
                  targets.use = targets.use, signaling = signaling, 
                  pairLR.use = pairLR.use, thresh = thresh, datasets = datasets, 
                  ligand.pvalues = ligand.pvalues, ligand.logFC = ligand.logFC, 
                  ligand.pct.1 = ligand.pct.1, ligand.pct.2 = ligand.pct.2, 
                  receptor.pvalues = receptor.pvalues, receptor.logFC = receptor.logFC, 
                  receptor.pct.1 = receptor.pct.1, receptor.pct.2 = receptor.pct.2, ...)
            }
        }
        else {
            LR <- data.frame()
            for (i in 1:length(object@LR)) {
                LR <- rbind(LR, object@LR[[i]]$LRsig)
            }
            LR <- unique(LR)
            cells.level <- levels(object@idents$joint)
            df.net <- subsetCommunication_internal(net, LR, cells.level, 
                slot.name = slot.name, sources.use = sources.use, 
                targets.use = targets.use, signaling = signaling, 
                pairLR.use = pairLR.use, thresh = thresh, datasets = datasets, 
                ligand.pvalues = ligand.pvalues, ligand.logFC = ligand.logFC, 
                ligand.pct.1 = ligand.pct.1, ligand.pct.2 = ligand.pct.2, 
                receptor.pvalues = receptor.pvalues, receptor.logFC = receptor.logFC, 
                receptor.pct.1 = receptor.pct.1, receptor.pct.2 = receptor.pct.2, ...)
        }
    }
    return(df.net)
}
environment(subsetCommunication) <- asNamespace('CellChat')
assignInNamespace("subsetCommunication", subsetCommunication, ns = "CellChat")
#rm(subsetCommunication)

netVisual_diffInteraction <- function (object, comparison = c(1, 2), measure = c("count", 
    "weight", "count.merged", "weight.merged"), color.use = NULL, 
    color.edge = c("#b2182b", "#2166ac"), title.name = NULL, 
    sources.use = NULL, targets.use = NULL, remove.isolate = FALSE, 
    top = 1, weight.scale = FALSE, vertex.weight = 20, vertex.weight.max = NULL, 
    vertex.size.max = 15, vertex.label.cex = 1, vertex.label.color = "black", 
    edge.weight.max = NULL, edge.width.max = 8, alpha.edge = 0.6, 
    label.edge = FALSE, edge.label.color = "black", edge.label.cex = 0.8, 
    edge.curved = 0.2, shape = "circle", layout = in_circle(), 
    margin = 0.2, arrow.width = 1, arrow.size = 0.2, net.diff = NULL) 
{
    options(warn = -1)
    measure <- match.arg(measure)
	if ( is.null(net.diff) ) { # here
    obj1 <- object@net[[comparison[1]]][[measure]]
    obj2 <- object@net[[comparison[2]]][[measure]]
    net.diff <- obj2 - obj1
	}
    if (measure %in% c("count", "count.merged")) {
        if (is.null(title.name)) {
            title.name = "Differential number of interactions"
        }
    }
    else if (measure %in% c("weight", "weight.merged")) {
        if (is.null(title.name)) {
            title.name = "Differential interaction strength"
        }
    }
    net <- net.diff
    if ((!is.null(sources.use)) | (!is.null(targets.use))) {
        df.net <- reshape2::melt(net, value.name = "value")
        colnames(df.net)[1:2] <- c("source", "target")
        if (!is.null(sources.use)) {
            if (is.numeric(sources.use)) {
                sources.use <- rownames(net.diff)[sources.use]
            }
            df.net <- subset(df.net, source %in% sources.use)
        }
        if (!is.null(targets.use)) {
            if (is.numeric(targets.use)) {
                targets.use <- rownames(net.diff)[targets.use]
            }
            df.net <- subset(df.net, target %in% targets.use)
        }
        cells.level <- rownames(net.diff)
        df.net$source <- factor(df.net$source, levels = cells.level)
        df.net$target <- factor(df.net$target, levels = cells.level)
        df.net$value[is.na(df.net$value)] <- 0
        net <- tapply(df.net[["value"]], list(df.net[["source"]], 
            df.net[["target"]]), sum)
    }
    if (remove.isolate) {
        idx1 <- which(Matrix::rowSums(net) == 0)
        idx2 <- which(Matrix::colSums(net) == 0)
        idx <- intersect(idx1, idx2)
        net <- net[-idx, ]
        net <- net[, -idx]
    }
    net[abs(net) < stats::quantile(abs(net), probs = 1 - top)] <- 0
    g <- graph_from_adjacency_matrix(net, mode = "directed", 
        weighted = T)
    edge.start <- igraph::ends(g, es = igraph::E(g), names = FALSE)
    coords <- layout_(g, layout)
    if (nrow(coords) != 1) {
        coords_scale = scale(coords)
    }
    else {
        coords_scale <- coords
    }
    if (is.null(color.use)) {
        color.use = scPalette(length(igraph::V(g)))
    }
    if (is.null(vertex.weight.max)) {
        vertex.weight.max <- max(vertex.weight)
    }
    vertex.weight <- vertex.weight/vertex.weight.max * vertex.size.max + 
        5
    loop.angle <- ifelse(coords_scale[igraph::V(g), 1] > 0, -atan(coords_scale[igraph::V(g), 
        2]/coords_scale[igraph::V(g), 1]), pi - atan(coords_scale[igraph::V(g), 
        2]/coords_scale[igraph::V(g), 1]))
    igraph::V(g)$size <- vertex.weight
    igraph::V(g)$color <- color.use[igraph::V(g)]
    igraph::V(g)$frame.color <- color.use[igraph::V(g)]
    igraph::V(g)$label.color <- vertex.label.color
    igraph::V(g)$label.cex <- vertex.label.cex
    if (label.edge) {
        igraph::E(g)$label <- igraph::E(g)$weight
        igraph::E(g)$label <- round(igraph::E(g)$label, digits = 1)
    }
    igraph::E(g)$arrow.width <- arrow.width
    igraph::E(g)$arrow.size <- arrow.size
    igraph::E(g)$label.color <- edge.label.color
    igraph::E(g)$label.cex <- edge.label.cex
    igraph::E(g)$color <- ifelse(igraph::E(g)$weight > 0, color.edge[1], 
        color.edge[2])
    igraph::E(g)$color <- grDevices::adjustcolor(igraph::E(g)$color, 
        alpha.edge)
	if ( ! is.null(igraph::E(g)$weight) ) # here
    igraph::E(g)$weight <- abs(igraph::E(g)$weight)
    if (is.null(edge.weight.max)) {
        edge.weight.max <- max(igraph::E(g)$weight)
    }
    if (weight.scale == TRUE) {
        igraph::E(g)$width <- 0.3 + igraph::E(g)$weight/edge.weight.max * 
            edge.width.max
    }
    else {
        igraph::E(g)$width <- 0.3 + edge.width.max * igraph::E(g)$weight
    }
    if (sum(edge.start[, 2] == edge.start[, 1]) != 0) {
        igraph::E(g)$loop.angle[which(edge.start[, 2] == edge.start[, 
            1])] <- loop.angle[edge.start[which(edge.start[, 
            2] == edge.start[, 1]), 1]]
    }
    radian.rescale <- function(x, start = 0, direction = 1) {
        c.rotate <- function(x) (x + start)%%(2 * pi) * direction
        c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
    }
    label.locs <- radian.rescale(x = 1:length(igraph::V(g)), 
        direction = -1, start = 0)
    label.dist <- vertex.weight/max(vertex.weight) + 2
    plot(g, edge.curved = edge.curved, vertex.shape = shape, 
        layout = coords_scale, margin = margin, vertex.label.dist = label.dist, 
        vertex.label.degree = label.locs, vertex.label.family = "Helvetica", 
        edge.label.family = "Helvetica")
    if (!is.null(title.name)) {
        text(0, 1.5, title.name, cex = 1.1)
    }
    gg <- recordPlot()
    return(gg)
}
environment(netVisual_diffInteraction) <- asNamespace('CellChat')
assignInNamespace("netVisual_diffInteraction", netVisual_diffInteraction, ns = "CellChat")

netVisual_heatmap <- function (object, comparison = c(1, 2), measure = c("count", 
    "weight"), signaling = NULL, slot.name = c("netP", "net"), 
    color.use = NULL, color.heatmap = c("#2166ac", "#b2182b"), 
    title.name = NULL, width = NULL, height = NULL, font.size = 8, 
    font.size.title = 10, cluster.rows = FALSE, cluster.cols = FALSE, 
    sources.use = NULL, targets.use = NULL, remove.isolate = FALSE, 
    row.show = NULL, col.show = NULL, legend.name = NULL) 
{
    if (!is.null(measure)) {
        measure <- match.arg(measure)
    }
    slot.name <- match.arg(slot.name)
    if (is.list(object@net[[1]])) {
        message("Do heatmap based on a merged object \n")
        obj1 <- object@net[[comparison[1]]][[measure]]
        obj2 <- object@net[[comparison[2]]][[measure]]
        net.diff <- obj2 - obj1
        if (measure == "count") {
            if (is.null(title.name)) {
                title.name = "Differential number of interactions"
            }
        }
        else if (measure == "weight") {
            if (is.null(title.name)) {
                title.name = "Differential interaction strength"
            }
        }
		if ( is.null(legend.name) ) # here
	        legend.name = "Relative values"
    }
    else {
        message("Do heatmap based on a single object \n")
        if (!is.null(signaling)) {
            net.diff <- slot(object, slot.name)$prob[, , signaling]
            if (is.null(title.name)) {
                title.name = paste0(signaling, " signaling network")
            }
			if ( is.null(legend.name) ) # here
	            legend.name <- "Communication Prob."
        }
        else if (!is.null(measure)) {
            net.diff <- object@net[[measure]]
            if (measure == "count") {
                if (is.null(title.name)) {
                  title.name = "Number of interactions"
                }
            }
            else if (measure == "weight") {
                if (is.null(title.name)) {
                  title.name = "Interaction strength"
                }
            }
			if ( is.null(legend.name) ) # here
	            legend.name <- title.name
        }
    }
    net <- net.diff
    if ((!is.null(sources.use)) | (!is.null(targets.use))) {
        df.net <- reshape2::melt(net, value.name = "value")
        colnames(df.net)[1:2] <- c("source", "target")
        if (!is.null(sources.use)) {
            if (is.numeric(sources.use)) {
                sources.use <- rownames(net.diff)[sources.use]
            }
            df.net <- subset(df.net, source %in% sources.use)
        }
        if (!is.null(targets.use)) {
            if (is.numeric(targets.use)) {
                targets.use <- rownames(net.diff)[targets.use]
            }
            df.net <- subset(df.net, target %in% targets.use)
        }
        cells.level <- rownames(net.diff)
        df.net$source <- factor(df.net$source, levels = cells.level)
        df.net$target <- factor(df.net$target, levels = cells.level)
        df.net$value[is.na(df.net$value)] <- 0
        net <- tapply(df.net[["value"]], list(df.net[["source"]], 
            df.net[["target"]]), sum)
    }
    net[is.na(net)] <- 0
    if (remove.isolate) {
        idx1 <- which(Matrix::rowSums(net) == 0)
        idx2 <- which(Matrix::colSums(net) == 0)
        idx <- intersect(idx1, idx2)
        net <- net[-idx, ]
        net <- net[, -idx]
        if (length(idx) > 0) {
            net <- net[-idx, ]
            net <- net[, -idx]
        }
    }
    mat <- net
    if (is.null(color.use)) {
        color.use <- scPalette(ncol(mat))
    }
    names(color.use) <- colnames(mat)
    if (!is.null(row.show)) {
        mat <- mat[row.show, ]
    }
    if (!is.null(col.show)) {
        mat <- mat[, col.show]
        color.use <- color.use[col.show]
    }
    if (min(mat) < 0) {
        color.heatmap.use = colorRamp3(c(min(mat), 0, max(mat)), 
            c(color.heatmap[1], "#f7f7f7", color.heatmap[2]))
        colorbar.break <- c(round(min(mat, na.rm = T), digits = nchar(sub(".*\\.(0*).*", 
            "\\1", min(mat, na.rm = T))) + 1), 0, round(max(mat, 
            na.rm = T), digits = nchar(sub(".*\\.(0*).*", "\\1", 
            max(mat, na.rm = T))) + 1))
    }
    else {
        if (length(color.heatmap) == 3) {
            color.heatmap.use = colorRamp3(c(0, min(mat), max(mat)), 
                color.heatmap)
        }
        else if (length(color.heatmap) == 2) {
            color.heatmap.use = colorRamp3(c(min(mat), max(mat)), 
                color.heatmap)
        }
        else if (length(color.heatmap) == 1) {
            color.heatmap.use = (grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, 
                name = color.heatmap))))(100)
        }
        colorbar.break <- c(round(min(mat, na.rm = T), digits = nchar(sub(".*\\.(0*).*", 
            "\\1", min(mat, na.rm = T))) + 1), round(max(mat, 
            na.rm = T), digits = nchar(sub(".*\\.(0*).*", "\\1", 
            max(mat, na.rm = T))) + 1))
    }
    df <- data.frame(group = colnames(mat))
    rownames(df) <- colnames(mat)
    col_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use), 
        which = "column", show_legend = FALSE, show_annotation_name = FALSE, 
        simple_anno_size = grid::unit(0.2, "cm"))
    row_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use), 
        which = "row", show_legend = FALSE, show_annotation_name = FALSE, 
        simple_anno_size = grid::unit(0.2, "cm"))
    ha1 = rowAnnotation(Strength = anno_barplot(rowSums(abs(mat)), 
        border = FALSE, gp = gpar(fill = color.use, col = color.use)), 
        show_annotation_name = FALSE)
    ha2 = HeatmapAnnotation(Strength = anno_barplot(colSums(abs(mat)), 
        border = FALSE, gp = gpar(fill = color.use, col = color.use)), 
        show_annotation_name = FALSE)
    if (sum(abs(mat) > 0) == 1) {
#        color.heatmap.use = c("white", color.heatmap.use) # here
    }
    else {
        mat[mat == 0] <- NA
    }
    ht1 = Heatmap(mat, col = color.heatmap.use, na_col = "white", 
        name = legend.name, bottom_annotation = col_annotation, 
        left_annotation = row_annotation, top_annotation = ha2, 
        right_annotation = ha1, cluster_rows = cluster.rows, 
        cluster_columns = cluster.rows, row_names_side = "left", 
        row_names_rot = 0, row_names_gp = gpar(fontsize = font.size), 
        column_names_gp = gpar(fontsize = font.size),
		column_title = 'Targets', column_title_side = 'bottom', # here 
        column_title_gp = gpar(fontsize = font.size.title), column_names_rot = 90, 
        row_title = "Sources", # here 
		row_title_gp = gpar(fontsize = font.size.title), 
        row_title_rot = 90, heatmap_legend_param = list(title_gp = gpar(fontsize = 8, 
            fontface = "plain"), title_position = "leftcenter-rot", 
            border = NA, legend_height = unit(20, "mm"), labels_gp = gpar(fontsize = 8), 
            grid_width = unit(2, "mm")))
#	ht1 = draw(ht1, column_title = title.name) # here
    return(ht1)
}
environment(netVisual_heatmap) <- asNamespace('CellChat')
assignInNamespace("netVisual_heatmap", netVisual_heatmap, ns = "CellChat")

netAnalysis_signalingChanges_scatter <- function (object, idents.use, color.use = c("grey10", "#F8766D", 
    "#00BFC4"), comparison = c(1, 2), signaling = NULL, signaling.label = NULL, 
    top.label = 1, signaling.exclude = NULL, xlims = NULL, ylims = NULL, 
    slot.name = "netP", dot.size = 2.5, point.shape = c(21, 22, 
        24, 23), label.size = 3, dot.alpha = 0.6, x.measure = "outdeg", 
    y.measure = "indeg", xlabel = "Differential outgoing interaction strength", 
    ylabel = "Differential incoming interaction strength", title = NULL, 
    font.size = 10, font.size.title = 10, do.label = T, show.legend = T, 
    show.axes = T) 
{
    if (is.list(object)) {
        object <- mergeCellChat(object, add.names = names(object))
    }
    if (is.list(object@net[[1]])) {
        dataset.name <- names(object@net)
        message(paste0("Visualizing differential outgoing and incoming signaling changes from ", 
            dataset.name[comparison[1]], " to ", dataset.name[comparison[2]]))
        title <- paste0("Signaling changes of ", idents.use, 
            " (", dataset.name[comparison[1]], " vs. ", dataset.name[comparison[2]], 
            ")")
        cell.levels <- levels(object@idents$joint)
        if (is.null(xlabel) | is.null(ylabel)) {
            xlabel = "Differential outgoing interaction strength"
            ylabel = "Differential incoming interaction strength"
        }
    }
    else {
        message("Visualizing outgoing and incoming signaling on a single object \n")
        title <- paste0("Signaling patterns of ", idents.use)
        if (length(slot(object, slot.name)$centr) == 0) {
            stop("Please run `netAnalysis_computeCentrality` to compute the network centrality scores! ")
        }
        cell.levels <- levels(object@idents)
    }
    if (!(idents.use %in% cell.levels)) {
        stop("Please check the input cell group names!")
    }
    if (is.null(signaling)) {
        signaling <- union(object@netP[[comparison[1]]]$pathways, 
            object@netP[[comparison[2]]]$pathways)
    }
    if (!is.null(signaling.exclude)) {
        signaling <- setdiff(signaling, signaling.exclude)
    }
    mat.all.merged <- list()
    for (ii in 1:length(comparison)) {
        if (length(slot(object, slot.name)[[comparison[ii]]]$centr) == 
            0) {
            stop("Please run `netAnalysis_computeCentrality` to compute the network centrality scores for each dataset seperately! ")
        }
        if (sum(c(x.measure, y.measure) %in% names(slot(object, 
            slot.name)[[comparison[ii]]]$centr[[1]])) != 2) {
            stop(paste0("`x.measure, y.measure` should be one of ", 
                paste(names(slot(object, slot.name)[[comparison[ii]]]$centr[[1]]), 
                  collapse = ", "), "\n", "`outdeg_unweighted` is only supported for version >= 1.1.2"))
        }
        centr <- slot(object, slot.name)[[comparison[ii]]]$centr
        outgoing <- matrix(0, nrow = length(cell.levels), ncol = length(centr))
        incoming <- matrix(0, nrow = length(cell.levels), ncol = length(centr))
        dimnames(outgoing) <- list(cell.levels, names(centr))
        dimnames(incoming) <- dimnames(outgoing)
        for (i in 1:length(centr)) {
            outgoing[, i] <- centr[[i]][[x.measure]]
            incoming[, i] <- centr[[i]][[y.measure]]
        }
        mat.out <- t(outgoing)
        mat.in <- t(incoming)
        mat.all <- array(0, dim = c(length(signaling), ncol(mat.out), 
            2))
        mat.t <- list(mat.out, mat.in)
        for (i in 1:length(comparison)) {
            mat = mat.t[[i]]
            mat1 <- mat[rownames(mat) %in% signaling, , drop = FALSE]
            mat <- matrix(0, nrow = length(signaling), ncol = ncol(mat))
            idx <- match(rownames(mat1), signaling)
            mat[idx[!is.na(idx)], ] <- mat1
            dimnames(mat) <- list(signaling, colnames(mat1))
            mat.all[, , i] = mat
        }
        dimnames(mat.all) <- list(dimnames(mat)[[1]], dimnames(mat)[[2]], 
            c("outgoing", "incoming"))
        mat.all.merged[[ii]] <- mat.all
    }
    mat.all.merged.use <- list(mat.all.merged[[1]][, idents.use, 
        ], mat.all.merged[[2]][, idents.use, ])
    idx.specific <- mat.all.merged.use[[1]] * mat.all.merged.use[[2]]
    mat.sum <- mat.all.merged.use[[2]] + mat.all.merged.use[[1]]
    out.specific.signaling <- rownames(idx.specific)[(mat.sum[, 
        1] != 0) & (idx.specific[, 1] == 0)]
    in.specific.signaling <- rownames(idx.specific)[(mat.sum[, 
        2] != 0) & (idx.specific[, 2] == 0)]
    mat.diff <- mat.all.merged.use[[2]] - mat.all.merged.use[[1]]
    idx <- rowSums(mat.diff) != 0
    mat.diff <- mat.diff[idx, , drop = FALSE] # here
    out.specific.signaling <- rownames(mat.diff) %in% out.specific.signaling
    in.specific.signaling <- rownames(mat.diff) %in% in.specific.signaling
    out.in.specific.signaling <- as.logical(out.specific.signaling * 
        in.specific.signaling)
    specificity.out.in <- matrix(0, nrow = nrow(mat.diff), ncol = 1)
    specificity.out.in[out.in.specific.signaling] <- 2
    specificity.out.in[setdiff(which(out.specific.signaling), 
        which(out.in.specific.signaling))] <- 1
    specificity.out.in[setdiff(which(in.specific.signaling), 
        which(out.in.specific.signaling))] <- -1
    df <- as.data.frame(mat.diff)
    df$specificity.out.in <- specificity.out.in
    df$specificity = matrix(0, nrow = nrow(mat.diff), ncol = 1) # here
    df$specificity[(specificity.out.in != 0) & (rowSums(mat.diff >= 
        0) == 2)] = 1
    df$specificity[(specificity.out.in != 0) & (rowSums(mat.diff <= 
        0) == 2)] = -1
    out.in.category <- c("Shared", "Incoming specific", "Outgoing specific", 
        "Incoming & Outgoing specific")
    specificity.category <- c("Shared", paste0(dataset.name[comparison[1]], 
        " specific"), paste0(dataset.name[comparison[2]], " specific"))
    df$specificity.out.in <- plyr::mapvalues(df$specificity.out.in, 
        from = c(0, -1, 1, 2), to = out.in.category)
    df$specificity.out.in <- factor(df$specificity.out.in, levels = out.in.category)
    df$specificity <- plyr::mapvalues(df$specificity, from = c(0, 
        -1, 1), to = specificity.category)
    df$specificity <- factor(df$specificity, levels = specificity.category)
    point.shape.use <- point.shape[out.in.category %in% unique(df$specificity.out.in)]
    df$specificity.out.in = droplevels(df$specificity.out.in, 
        exclude = setdiff(out.in.category, unique(df$specificity.out.in)))
    color.use <- color.use[specificity.category %in% unique(df$specificity)]
    df$specificity = droplevels(df$specificity, exclude = setdiff(specificity.category, 
        unique(df$specificity)))
    df$labels <- rownames(df)
    gg <- ggplot(data = df, aes(outgoing, incoming)) + geom_point(aes(colour = specificity, 
        fill = specificity, shape = specificity.out.in), size = dot.size)
    gg <- gg + theme_linedraw() + theme(panel.grid = element_blank()) + 
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", 
            size = 0.25) + geom_vline(xintercept = 0, linetype = "dashed", 
        color = "grey50", size = 0.25) + theme(text = element_text(size = font.size), 
        legend.key.height = grid::unit(0.15, "in")) + labs(title = title, 
        x = xlabel, y = ylabel) + theme(plot.title = element_text(size = font.size.title, 
        hjust = 0.5, face = "plain")) + theme(axis.line.x = element_line(size = 0.25), 
        axis.line.y = element_line(size = 0.25))
    gg <- gg + scale_fill_manual(values = ggplot2::alpha(color.use, 
        alpha = dot.alpha), drop = FALSE) + guides(fill = "none")
    gg <- gg + scale_colour_manual(values = color.use, drop = FALSE)
    gg <- gg + scale_shape_manual(values = point.shape.use)
    gg <- gg + theme(legend.title = element_blank())
    if (!is.null(xlims)) {
        gg <- gg + xlim(xlims)
    }
    if (!is.null(ylims)) {
        gg <- gg + ylim(ylims)
    }
    if (do.label) {
        if (is.null(signaling.label)) {
            thresh <- stats::quantile(abs(as.matrix(df[, 1:2])), 
                probs = 1 - top.label)
            idx = abs(df[, 1]) > thresh | abs(df[, 2]) > thresh
            data.label <- df[idx, ]
        }
        else {
            data.label <- df[rownames(df) %in% signaling.label, 
                ]
        }
        gg <- gg + ggrepel::geom_text_repel(data = data.label, 
            mapping = aes(label = labels, colour = specificity), 
            size = label.size, show.legend = F, segment.size = 0.2, 
            segment.alpha = 0.5)
    }
    if (!show.legend) {
        gg <- gg + theme(legend.position = "none")
    }
    if (!show.axes) {
        gg <- gg + theme_void()
    }
    gg
}
environment(netAnalysis_signalingChanges_scatter) <- asNamespace('CellChat')
assignInNamespace("netAnalysis_signalingChanges_scatter", netAnalysis_signalingChanges_scatter, ns = "CellChat")



netVisual_bubble <- function (object, sources.use = NULL, targets.use = NULL, signaling = NULL, 
    pairLR.use = NULL, color.heatmap = c("Spectral", "viridis"), 
    n.colors = 10, direction = -1, thresh = 0.05, comparison = NULL, 
    group = NULL, remove.isolate = FALSE, max.dataset = NULL, 
    min.dataset = NULL, min.quantile = 0, max.quantile = 1, line.on = TRUE, 
    line.size = 0.2, color.text.use = TRUE, color.text = NULL, 
    title.name = NULL, font.size = 10, font.size.title = 10, 
    show.legend = TRUE, grid.on = TRUE, color.grid = "grey90", 
    angle.x = 90, vjust.x = NULL, hjust.x = NULL, return.data = FALSE) 
{
    color.heatmap <- match.arg(color.heatmap)
    if (is.list(object@net[[1]])) {
        message("Comparing communications on a merged object \n")
    }
    else {
        message("Comparing communications on a single object \n")
    }
    if (is.null(vjust.x) | is.null(hjust.x)) {
        angle = c(0, 45, 90)
        hjust = c(0, 1, 1)
        vjust = c(0, 1, 0.5)
        vjust.x = vjust[angle == angle.x]
        hjust.x = hjust[angle == angle.x]
    }
    if (length(color.heatmap) == 1) {
        color.use <- tryCatch({
            RColorBrewer::brewer.pal(n = n.colors, name = color.heatmap)
        }, error = function(e) {
            (scales::viridis_pal(option = color.heatmap, direction = -1))(n.colors)
        })
    }
    else {
        color.use <- color.heatmap
    }
    if (direction == -1) {
        color.use <- rev(color.use)
    }
    if (is.null(comparison)) {
        cells.level <- levels(object@idents)
        if (is.numeric(sources.use)) {
            sources.use <- cells.level[sources.use]
        }
        if (is.numeric(targets.use)) {
            targets.use <- cells.level[targets.use]
        }
        df.net <- subsetCommunication(object, slot.name = "net", 
            sources.use = sources.use, targets.use = targets.use, 
            signaling = signaling, pairLR.use = pairLR.use, thresh = thresh)
        df.net$source.target <- paste(df.net$source, df.net$target, 
            sep = " -> ")
        source.target <- paste(rep(sources.use, each = length(targets.use)), 
            targets.use, sep = " -> ")
        source.target.isolate <- setdiff(source.target, unique(df.net$source.target))
        if (length(source.target.isolate) > 0) {
            df.net.isolate <- as.data.frame(matrix(NA, nrow = length(source.target.isolate), 
                ncol = ncol(df.net)))
            colnames(df.net.isolate) <- colnames(df.net)
            df.net.isolate$source.target <- source.target.isolate
            df.net.isolate$interaction_name_2 <- df.net$interaction_name_2[1]
            df.net.isolate$pval <- 1
            a <- stringr::str_split(df.net.isolate$source.target, 
                " -> ", simplify = T)
            df.net.isolate$source <- as.character(a[, 1])
            df.net.isolate$target <- as.character(a[, 2])
            df.net <- rbind(df.net, df.net.isolate)
        }
		df.net$pval.original <- df.net$pval # here
        df.net$pval[df.net$pval > 0.05] = 1
        df.net$pval[df.net$pval > 0.01 & df.net$pval <= 0.05] = 2
        df.net$pval[df.net$pval <= 0.01] = 3
        df.net$prob[df.net$prob == 0] <- NA
        df.net$prob.original <- df.net$prob
        df.net$prob <- -1/log(df.net$prob)
        idx1 <- which(is.infinite(df.net$prob) | df.net$prob < 
            0)
        if (sum(idx1) > 0) {
            values.assign <- seq(max(df.net$prob, na.rm = T) * 
                1.1, max(df.net$prob, na.rm = T) * 1.5, length.out = length(idx1))
            position <- sort(prob.original[idx1], index.return = TRUE)$ix
            df.net$prob[idx1] <- values.assign[match(1:length(idx1), 
                position)]
        }
        df.net$source <- factor(df.net$source, levels = cells.level[cells.level %in% 
            unique(df.net$source)])
        df.net$target <- factor(df.net$target, levels = cells.level[cells.level %in% 
            unique(df.net$target)])
        group.names <- paste(rep(levels(df.net$source), each = length(levels(df.net$target))), 
            levels(df.net$target), sep = " -> ")
        df.net$interaction_name_2 <- as.character(df.net$interaction_name_2)
        df.net <- with(df.net, df.net[order(interaction_name_2), 
            ])
        df.net$interaction_name_2 <- factor(df.net$interaction_name_2, 
            levels = unique(df.net$interaction_name_2))
        cells.order <- group.names
        df.net$source.target <- factor(df.net$source.target, 
            levels = cells.order)
        df <- df.net
    }
    else {
        dataset.name <- names(object@net)
        df.net.all <- subsetCommunication(object, slot.name = "net", 
            sources.use = sources.use, targets.use = targets.use, 
            signaling = signaling, pairLR.use = pairLR.use, thresh = thresh)
        df.all <- data.frame()
        for (ii in 1:length(comparison)) {
            cells.level <- levels(object@idents[[comparison[ii]]])
            if (is.numeric(sources.use)) {
                sources.use <- cells.level[sources.use]
            }
            if (is.numeric(targets.use)) {
                targets.use <- cells.level[targets.use]
            }
            df.net <- df.net.all[[comparison[ii]]]
            df.net$interaction_name_2 <- as.character(df.net$interaction_name_2)
            df.net$source.target <- paste(df.net$source, df.net$target, 
                sep = " -> ")
            source.target <- paste(rep(sources.use, each = length(targets.use)), 
                targets.use, sep = " -> ")
            source.target.isolate <- setdiff(source.target, unique(df.net$source.target))
            if (length(source.target.isolate) > 0) {
                df.net.isolate <- as.data.frame(matrix(NA, nrow = length(source.target.isolate), 
                  ncol = ncol(df.net)))
                colnames(df.net.isolate) <- colnames(df.net)
                df.net.isolate$source.target <- source.target.isolate
                df.net.isolate$interaction_name_2 <- df.net$interaction_name_2[1]
                df.net.isolate$pval <- 1
                a <- stringr::str_split(df.net.isolate$source.target, 
                  " -> ", simplify = T)
                df.net.isolate$source <- as.character(a[, 1])
                df.net.isolate$target <- as.character(a[, 2])
                df.net <- rbind(df.net, df.net.isolate)
            }
            df.net$source <- factor(df.net$source, levels = cells.level[cells.level %in% 
                unique(df.net$source)])
            df.net$target <- factor(df.net$target, levels = cells.level[cells.level %in% 
                unique(df.net$target)])
            group.names <- paste(rep(levels(df.net$source), each = length(levels(df.net$target))), 
                levels(df.net$target), sep = " -> ")
            group.names0 <- group.names
            group.names <- paste0(group.names0, " (", dataset.name[comparison[ii]], 
                ")")
            if (nrow(df.net) > 0) {
				df.net$pval.original <- df.net$pval ## here
                df.net$pval[df.net$pval > 0.05] = 1
                df.net$pval[df.net$pval > 0.01 & df.net$pval <= 
                  0.05] = 2
                df.net$pval[df.net$pval <= 0.01] = 3
                df.net$prob[df.net$prob == 0] <- NA
                df.net$prob.original <- df.net$prob
                df.net$prob <- -1/log(df.net$prob)
            }
            else {
                df.net <- as.data.frame(matrix(NA, nrow = length(group.names), 
                  ncol = 5))
                colnames(df.net) <- c("interaction_name_2", "source.target", 
                  "prob", "pval", "prob.original")
                df.net$source.target <- group.names0
            }
            df.net$group.names <- as.character(df.net$source.target)
            df.net$source.target <- paste0(df.net$source.target, 
                " (", dataset.name[comparison[ii]], ")")
            df.net$dataset <- dataset.name[comparison[ii]]
            df.all <- rbind(df.all, df.net)
        }
        if (nrow(df.all) == 0) {
            stop("No interactions are detected. Please consider changing the cell groups for analysis. ")
        }
        idx1 <- which(is.infinite(df.all$prob) | df.all$prob < 
            0)
        if (sum(idx1) > 0) {
            values.assign <- seq(max(df.all$prob, na.rm = T) * 
                1.1, max(df.all$prob, na.rm = T) * 1.5, length.out = length(idx1))
            position <- sort(df.all$prob.original[idx1], index.return = TRUE)$ix
            df.all$prob[idx1] <- values.assign[match(1:length(idx1), 
                position)]
        }
        df.all$interaction_name_2[is.na(df.all$interaction_name_2)] <- df.all$interaction_name_2[!is.na(df.all$interaction_name_2)][1]
        df <- df.all
        df <- with(df, df[order(interaction_name_2), ])
        df$interaction_name_2 <- factor(df$interaction_name_2, 
            levels = unique(df$interaction_name_2))
        cells.order <- c()
        dataset.name.order <- c()
        for (i in 1:length(group.names0)) {
            for (j in 1:length(comparison)) {
                cells.order <- c(cells.order, paste0(group.names0[i], 
                  " (", dataset.name[comparison[j]], ")"))
                dataset.name.order <- c(dataset.name.order, dataset.name[comparison[j]])
            }
        }
        df$source.target <- factor(df$source.target, levels = cells.order)
    }
    min.cutoff <- quantile(df$prob, min.quantile, na.rm = T)
    max.cutoff <- quantile(df$prob, max.quantile, na.rm = T)
    df$prob[df$prob < min.cutoff] <- min.cutoff
    df$prob[df$prob > max.cutoff] <- max.cutoff
    if (remove.isolate) {
        df <- df[!is.na(df$prob), ]
        line.on <- FALSE
    }
    if (!is.null(max.dataset)) {
        signaling <- as.character(unique(df$interaction_name_2))
        for (i in signaling) {
            df.i <- df[df$interaction_name_2 == i, , drop = FALSE]
            cell <- as.character(unique(df.i$group.names))
            for (j in cell) {
                df.i.j <- df.i[df.i$group.names == j, , drop = FALSE]
                values <- df.i.j$prob
                idx.max <- which(values == max(values, na.rm = T))
                idx.min <- which(values == min(values, na.rm = T))
                dataset.na <- c(df.i.j$dataset[is.na(values)], 
                  setdiff(dataset.name[comparison], df.i.j$dataset))
				if (length(dataset.na) > 0 ) {
						if ( length(dataset.na) != 1 ) {
								stop()
						}
						if ( dataset.name[max.dataset] == dataset.na ) {
								df.i.j$prob <- NA
						}
				} else {
						if ( idx.max != max.dataset ) {
								df.i.j$prob <- NA
						}
				}
#                if (length(idx.max) > 0) {
#                  if (!(df.i.j$dataset[idx.max] %in% dataset.name[max.dataset])) {
#                    df.i.j$prob <- NA
#                  }
#                  else if ((idx.max != idx.min) & !is.null(min.dataset)) {
#                    if (!(df.i.j$dataset[idx.min] %in% dataset.name[min.dataset])) {
#                      df.i.j$prob <- NA
#                    }
#                    else if (length(dataset.na) > 0 & sum(!(dataset.name[min.dataset] %in% 
#                      dataset.na)) > 0) {
#                      df.i.j$prob <- NA
#                    }
#                  }
#                }
                df.i[df.i$group.names == j, "prob"] <- df.i.j$prob
            }
            df[df$interaction_name_2 == i, "prob"] <- df.i$prob
        }
    }
#	df2 <- df
#	df[df$prob == 'NA'] <- NA
    if (remove.isolate) {
        df.test <- df[!is.na(df$prob), ]
		if ( nrow(df.test) != 0 ) {
			df <- df.test
		}
        line.on <- FALSE
    }
    if (nrow(df) == 0) {
        stop("No interactions are detected. Please consider changing the cell groups for analysis. ")
    }
    df$interaction_name_2 <- factor(df$interaction_name_2, levels = unique(df$interaction_name_2))
    df$source.target = droplevels(df$source.target, exclude = setdiff(levels(df$source.target), 
        unique(df$source.target)))
    g <- ggplot(df, aes(x = source.target, y = interaction_name_2, 
        color = prob, size = pval)) + geom_point(pch = 16) + 
        theme_linedraw() + theme(panel.grid.major = element_blank()) + 
        theme(axis.text.x = element_text(angle = angle.x, hjust = hjust.x, 
            vjust = vjust.x), axis.title.x = element_blank(), 
            axis.title.y = element_blank()) + scale_x_discrete(position = "bottom")
    values <- c(1, 2, 3)
    names(values) <- c("p > 0.05", "0.01 < p < 0.05", "p < 0.01")
    g <- g + scale_radius(range = c(min(df$pval), max(df$pval)), 
        breaks = sort(unique(df$pval)), labels = names(values)[values %in% 
            sort(unique(df$pval))], name = "p-value")
    if (min(df$prob, na.rm = T) != max(df$prob, na.rm = T)) {
        g <- g + scale_colour_gradientn(colors = colorRampPalette(color.use)(99), 
            na.value = "white", limits = c(quantile(df$prob, 
                0, na.rm = T), quantile(df$prob, 1, na.rm = T)), 
            breaks = c(quantile(df$prob, 0, na.rm = T), quantile(df$prob, 
                1, na.rm = T)), labels = c("min", "max")) + guides(color = guide_colourbar(barwidth = 0.5, 
            title = "Commun. Prob.", order = 1)) # here
    }
    else {
        g <- g + scale_colour_gradientn(colors = colorRampPalette(color.use)(99), 
            na.value = "white") + guides(color = guide_colourbar(barwidth = 0.5, 
            title = "Commun. Prob.", order = 1)) # here
    }
    g <- g + theme(text = element_text(size = font.size), plot.title = element_text(size = font.size.title)) + 
        theme(legend.title = element_text(size = 8), legend.text = element_text(size = 6))
    if (grid.on) {
        if (length(unique(df$source.target)) > 1) {
            g <- g + geom_vline(xintercept = seq(1.5, length(unique(df$source.target)) - 
                0.5, 1), lwd = 0.1, colour = color.grid)
        }
        if (length(unique(df$interaction_name_2)) > 1) {
            g <- g + geom_hline(yintercept = seq(1.5, length(unique(df$interaction_name_2)) - 
                0.5, 1), lwd = 0.1, colour = color.grid)
        }
    }
    if (!is.null(title.name)) {
        g <- g + ggtitle(title.name) + theme(plot.title = element_text(hjust = 0.5))
    }
    if (!is.null(comparison)) {
        if (line.on) {
            xintercept = seq(0.5 + length(dataset.name[comparison]), 
                length(group.names0) * length(dataset.name[comparison]), 
                by = length(dataset.name[comparison]))
            g <- g + geom_vline(xintercept = xintercept, linetype = "dashed", 
                color = "grey60", size = line.size)
        }
        if (color.text.use) {
            if (is.null(group)) {
                group <- 1:length(comparison)
                names(group) <- dataset.name[comparison]
            }
            if (is.null(color.text)) {
                color <- ggPalette(length(unique(group)))
            }
            else {
                color <- color.text
            }
            names(color) <- names(group[!duplicated(group)])
            color <- color[group]
            dataset.name.order <- levels(df$source.target)
            dataset.name.order <- stringr::str_match(dataset.name.order, 
                "\\(.*\\)")
            dataset.name.order <- stringr::str_sub(dataset.name.order, 
                2, stringr::str_length(dataset.name.order) - 
                  1)
            xtick.color <- color[dataset.name.order]
            g <- g + theme(axis.text.x = element_text(colour = xtick.color))
        }
    }
    if (!show.legend) {
        g <- g + theme(legend.position = "none")
    }
    if (return.data) {
        return(list(communication = df, gg.obj = g))
    }
    else {
        return(g)
    }
}
environment(netVisual_bubble) <- asNamespace('CellChat')
assignInNamespace("netVisual_bubble", netVisual_bubble, ns = "CellChat")



netVisual_hierarchy1 <- function (net, vertex.receiver, color.use = NULL, title.name = NULL, 
    sources.use = NULL, targets.use = NULL, remove.isolate = FALSE, 
    top = 1, weight.scale = FALSE, vertex.weight = 20, vertex.weight.max = NULL, 
    vertex.size.max = NULL, edge.weight.max = NULL, edge.width.max = 8, 
    alpha.edge = 0.6, label.dist = 2.8, space.v = 1.5, space.h = 1.6, 
    shape = NULL, label.edge = FALSE, edge.curved = 0, margin = 0.2, 
    vertex.label.cex = 0.6, vertex.label.color = "black", arrow.width = 1, 
    arrow.size = 0.2, edge.label.color = "black", edge.label.cex = 0.5, 
    vertex.size = NULL) 
{
    if (!is.null(vertex.size)) {
        warning("'vertex.size' is deprecated. Use `vertex.weight`")
    }
    if (is.null(vertex.size.max)) {
        if (length(unique(vertex.weight)) == 1) {
            vertex.size.max <- 5
        }
        else {
            vertex.size.max <- 15
        }
    }
    options(warn = -1)
    thresh <- stats::quantile(net, probs = 1 - top)
    net[net < thresh] <- 0
    cells.level <- rownames(net)
    if ((!is.null(sources.use)) | (!is.null(targets.use))) {
        df.net <- reshape2::melt(net, value.name = "value")
        colnames(df.net)[1:2] <- c("source", "target")
        if (!is.null(sources.use)) {
            if (is.numeric(sources.use)) {
                sources.use <- cells.level[sources.use]
            }
            df.net <- subset(df.net, source %in% sources.use)
        }
        if (!is.null(targets.use)) {
            if (is.numeric(targets.use)) {
                targets.use <- cells.level[targets.use]
            }
            df.net <- subset(df.net, target %in% targets.use)
        }
        df.net$source <- factor(df.net$source, levels = cells.level)
        df.net$target <- factor(df.net$target, levels = cells.level)
        df.net$value[is.na(df.net$value)] <- 0
        net <- tapply(df.net[["value"]], list(df.net[["source"]], 
            df.net[["target"]]), sum)
    }
    net[is.na(net)] <- 0
    if (remove.isolate) {
        idx1 <- which(Matrix::rowSums(net) == 0)
        idx2 <- which(Matrix::colSums(net) == 0)
        idx <- intersect(idx1, idx2)
        net <- net[-idx, ]
        net <- net[, -idx]
    }
    if (is.null(color.use)) {
        color.use <- scPalette(nrow(net))
    }
    if (is.null(vertex.weight.max)) {
        vertex.weight.max <- max(vertex.weight)
    }
    vertex.weight <- vertex.weight/vertex.weight.max * vertex.size.max + 
        6
    m <- length(vertex.receiver)
    net2 <- net
    reorder.row <- c(vertex.receiver, setdiff(1:nrow(net), vertex.receiver))
    net2 <- net2[reorder.row, vertex.receiver, drop = FALSE] ## here
    m1 <- nrow(net2)
    n1 <- ncol(net2)
    net3 <- rbind(cbind(matrix(0, m1, m1), net2), matrix(0, n1, 
        m1 + n1))
    row.names(net3) <- c(row.names(net)[vertex.receiver], row.names(net)[setdiff(1:m1, 
        vertex.receiver)], rep("", m))
    colnames(net3) <- row.names(net3)
    color.use3 <- c(color.use[vertex.receiver], color.use[setdiff(1:m1, 
        vertex.receiver)], rep("#FFFFFF", length(vertex.receiver)))
    color.use3.frame <- c(color.use[vertex.receiver], color.use[setdiff(1:m1, 
        vertex.receiver)], color.use[vertex.receiver])
    if (length(vertex.weight) != 1) {
        vertex.weight = c(vertex.weight[vertex.receiver], vertex.weight[setdiff(1:m1, 
            vertex.receiver)], vertex.weight[vertex.receiver])
    }
    if (is.null(shape)) {
        shape <- c(rep("circle", m), rep("circle", m1 - m), rep("circle", 
            m))
    }
    g <- graph_from_adjacency_matrix(net3, mode = "directed", 
        weighted = T)
    edge.start <- ends(g, es = E(g), names = FALSE)
    coords <- matrix(NA, nrow(net3), 2)
    coords[1:m, 1] <- 0
    coords[(m + 1):m1, 1] <- space.h
    coords[(m1 + 1):nrow(net3), 1] <- space.h/2
    coords[1:m, 2] <- seq(space.v, 0, by = -space.v/(m - 1))
    coords[(m + 1):m1, 2] <- seq(space.v, 0, by = -space.v/(m1 - 
        m - 1))
    coords[(m1 + 1):nrow(net3), 2] <- seq(space.v, 0, by = -space.v/(n1 - 
        1))
    coords[is.nan(coords)] <- 1.5 ## here
    coords_scale <- coords
    igraph::V(g)$size <- vertex.weight
    igraph::V(g)$color <- color.use3[igraph::V(g)]
    igraph::V(g)$frame.color <- color.use3.frame[igraph::V(g)]
    igraph::V(g)$label.color <- vertex.label.color
    igraph::V(g)$label.cex <- vertex.label.cex
    if (label.edge) {
        E(g)$label <- E(g)$weight
        igraph::E(g)$label <- round(igraph::E(g)$label, digits = 1)
    }
    if (is.null(edge.weight.max)) {
        edge.weight.max <- max(igraph::E(g)$weight)
    }
    if (weight.scale == TRUE) {
        E(g)$width <- 0.3 + E(g)$weight/edge.weight.max * edge.width.max
    }
    else {
        E(g)$width <- 0.3 + edge.width.max * E(g)$weight
    }
    E(g)$arrow.width <- arrow.width
    E(g)$arrow.size <- arrow.size
    E(g)$label.color <- edge.label.color
    E(g)$label.cex <- edge.label.cex
    E(g)$color <- adjustcolor(igraph::V(g)$color[edge.start[, 
        1]], alpha.edge)
    label.dist <- c(rep(space.h * label.dist, m), rep(space.h * 
        label.dist, m1 - m), rep(0, nrow(net3) - m1))
    label.locs <- c(rep(-pi, m), rep(0, m1 - m), rep(-pi, nrow(net3) - 
        m1))
    text.pos <- cbind(c(-space.h/1.5, space.h/22, space.h/1.5), 
        space.v - space.v/7)
    igraph::add.vertex.shape("fcircle", clip = igraph::igraph.shape.noclip, 
        plot = mycircle, parameters = list(vertex.frame.color = 1, 
            vertex.frame.width = 1))
    plot(g, edge.curved = edge.curved, layout = coords_scale, 
        margin = margin, rescale = T, vertex.shape = "fcircle", 
        vertex.frame.width = c(rep(1, m1), rep(2, nrow(net3) - 
            m1)), vertex.label.degree = label.locs, vertex.label.dist = label.dist, 
        vertex.label.family = "Helvetica")
    text(text.pos, c("Source", "Target", "Source"), cex = 0.8, 
        col = c("#c51b7d", "#c51b7d", "#2f6661"))
    arrow.pos1 <- c(-space.h/1.5, space.v - space.v/4, space.h/1e+05, 
        space.v - space.v/4)
    arrow.pos2 <- c(space.h/1.5, space.v - space.v/4, space.h/20, 
        space.v - space.v/4)
    shape::Arrows(arrow.pos1[1], arrow.pos1[2], arrow.pos1[3], 
        arrow.pos1[4], col = "#c51b7d", arr.lwd = 1e-04, arr.length = 0.2, 
        lwd = 0.8, arr.type = "triangle")
    shape::Arrows(arrow.pos2[1], arrow.pos2[2], arrow.pos2[3], 
        arrow.pos2[4], col = "#2f6661", arr.lwd = 1e-04, arr.length = 0.2, 
        lwd = 0.8, arr.type = "triangle")
    if (!is.null(title.name)) {
        title.pos = c(space.h/8, space.v)
        text(title.pos[1], title.pos[2], paste0(title.name, " signaling network"), 
            cex = 1)
    }
    gg <- recordPlot()
    return(gg)
}
environment(netVisual_hierarchy1) <- asNamespace('CellChat')
assignInNamespace("netVisual_hierarchy1", netVisual_hierarchy1, ns = "CellChat")



netVisual_hierarchy2 <- function (net, vertex.receiver, color.use = NULL, title.name = NULL, 
    sources.use = NULL, targets.use = NULL, remove.isolate = FALSE, 
    top = 1, weight.scale = FALSE, vertex.weight = 20, vertex.weight.max = NULL, 
    vertex.size.max = NULL, edge.weight.max = NULL, edge.width.max = 8, 
    alpha.edge = 0.6, label.dist = 2.8, space.v = 1.5, space.h = 1.6, 
    shape = NULL, label.edge = FALSE, edge.curved = 0, margin = 0.2, 
    vertex.label.cex = 0.6, vertex.label.color = "black", arrow.width = 1, 
    arrow.size = 0.2, edge.label.color = "black", edge.label.cex = 0.5, 
    vertex.size = NULL) 
{
    if (!is.null(vertex.size)) {
        warning("'vertex.size' is deprecated. Use `vertex.weight`")
    }
    if (is.null(vertex.size.max)) {
        if (length(unique(vertex.weight)) == 1) {
            vertex.size.max <- 5
        }
        else {
            vertex.size.max <- 15
        }
    }
    options(warn = -1)
    thresh <- stats::quantile(net, probs = 1 - top)
    net[net < thresh] <- 0
    if ((!is.null(sources.use)) | (!is.null(targets.use))) {
        df.net <- reshape2::melt(net, value.name = "value")
        colnames(df.net)[1:2] <- c("source", "target")
        if (!is.null(sources.use)) {
            if (is.numeric(sources.use)) {
                sources.use <- levels(object@idents)[sources.use]
            }
            df.net <- subset(df.net, source %in% sources.use)
        }
        if (!is.null(targets.use)) {
            if (is.numeric(targets.use)) {
                targets.use <- levels(object@idents)[targets.use]
            }
            df.net <- subset(df.net, target %in% targets.use)
        }
        cells.level <- levels(object@idents)
        df.net$source <- factor(df.net$source, levels = cells.level)
        df.net$target <- factor(df.net$target, levels = cells.level)
        df.net$value[is.na(df.net$value)] <- 0
        net <- tapply(df.net[["value"]], list(df.net[["source"]], 
            df.net[["target"]]), sum)
    }
    net[is.na(net)] <- 0
    if (remove.isolate) {
        idx1 <- which(Matrix::rowSums(net) == 0)
        idx2 <- which(Matrix::colSums(net) == 0)
        idx <- intersect(idx1, idx2)
        net <- net[-idx, ]
        net <- net[, -idx]
    }
    if (is.null(color.use)) {
        color.use <- scPalette(nrow(net))
    }
    if (is.null(vertex.weight.max)) {
        vertex.weight.max <- max(vertex.weight)
    }
    vertex.weight <- vertex.weight/vertex.weight.max * vertex.size.max + 
        6
    m <- length(vertex.receiver)
    m0 <- nrow(net) - length(vertex.receiver)
    net2 <- net
    reorder.row <- c(setdiff(1:nrow(net), vertex.receiver), vertex.receiver)
    net2 <- net2[reorder.row, vertex.receiver, drop = FALSE] ## here
    m1 <- nrow(net2)
    n1 <- ncol(net2)
    net3 <- rbind(cbind(matrix(0, m1, m1), net2), matrix(0, n1, 
        m1 + n1))
    row.names(net3) <- c(row.names(net)[setdiff(1:m1, vertex.receiver)], 
        row.names(net)[vertex.receiver], rep("", m))
    colnames(net3) <- row.names(net3)
    color.use3 <- c(color.use[setdiff(1:m1, vertex.receiver)], 
        color.use[vertex.receiver], rep("#FFFFFF", length(vertex.receiver)))
    color.use3.frame <- c(color.use[setdiff(1:m1, vertex.receiver)], 
        color.use[vertex.receiver], color.use[vertex.receiver])
    if (length(vertex.weight) != 1) {
        vertex.weight = c(vertex.weight[setdiff(1:m1, vertex.receiver)], 
            vertex.weight[vertex.receiver], vertex.weight[vertex.receiver])
    }
    if (is.null(shape)) {
        shape <- rep("circle", nrow(net3))
    }
    g <- graph_from_adjacency_matrix(net3, mode = "directed", 
        weighted = T)
    edge.start <- ends(g, es = igraph::E(g), names = FALSE)
    coords <- matrix(NA, nrow(net3), 2)
    coords[1:m0, 1] <- 0
    coords[(m0 + 1):m1, 1] <- space.h
    coords[(m1 + 1):nrow(net3), 1] <- space.h/2
    coords[1:m0, 2] <- seq(space.v, 0, by = -space.v/(m0 - 1))
    coords[(m0 + 1):m1, 2] <- seq(space.v, 0, by = -space.v/(m1 - 
        m0 - 1))
    coords[(m1 + 1):nrow(net3), 2] <- seq(space.v, 0, by = -space.v/(n1 - 
        1))
    coords[is.nan(coords)] <- 1.5 ## here
    coords_scale <- coords
    igraph::V(g)$size <- vertex.weight
    igraph::V(g)$color <- color.use3[igraph::V(g)]
    igraph::V(g)$frame.color <- color.use3.frame[igraph::V(g)]
    igraph::V(g)$label.color <- vertex.label.color
    igraph::V(g)$label.cex <- vertex.label.cex
    if (label.edge) {
        igraph::E(g)$label <- igraph::E(g)$weight
        igraph::E(g)$label <- round(igraph::E(g)$label, digits = 1)
    }
    if (is.null(edge.weight.max)) {
        edge.weight.max <- max(igraph::E(g)$weight)
    }
    if (weight.scale == TRUE) {
        igraph::E(g)$width <- 0.3 + igraph::E(g)$weight/edge.weight.max * 
            edge.width.max
    }
    else {
        igraph::E(g)$width <- 0.3 + edge.width.max * igraph::E(g)$weight
    }
    igraph::E(g)$arrow.width <- arrow.width
    igraph::E(g)$arrow.size <- arrow.size
    igraph::E(g)$label.color <- edge.label.color
    igraph::E(g)$label.cex <- edge.label.cex
    igraph::E(g)$color <- adjustcolor(igraph::V(g)$color[edge.start[, 
        1]], alpha.edge)
    label.dist <- c(rep(space.h * label.dist, m), rep(space.h * 
        label.dist, m1 - m), rep(0, nrow(net3) - m1))
    label.locs <- c(rep(-pi, m0), rep(0, m1 - m0), rep(-pi, nrow(net3) - 
        m1))
    text.pos <- cbind(c(-space.h/1.5, space.h/22, space.h/1.5), 
        space.v - space.v/7)
    igraph::add.vertex.shape("fcircle", clip = igraph::igraph.shape.noclip, 
        plot = mycircle, parameters = list(vertex.frame.color = 1, 
            vertex.frame.width = 1))
    plot(g, edge.curved = edge.curved, layout = coords_scale, 
        margin = margin, rescale = T, vertex.shape = "fcircle", 
        vertex.frame.width = c(rep(1, m1), rep(2, nrow(net3) - 
            m1)), vertex.label.degree = label.locs, vertex.label.dist = label.dist, 
        vertex.label.family = "Helvetica")
    text(text.pos, c("Source", "Target", "Source"), cex = 0.8, 
        col = c("#c51b7d", "#2f6661", "#2f6661"))
    arrow.pos1 <- c(-space.h/1.5, space.v - space.v/4, space.h/1e+05, 
        space.v - space.v/4)
    arrow.pos2 <- c(space.h/1.5, space.v - space.v/4, space.h/20, 
        space.v - space.v/4)
    shape::Arrows(arrow.pos1[1], arrow.pos1[2], arrow.pos1[3], 
        arrow.pos1[4], col = "#c51b7d", arr.lwd = 1e-04, arr.length = 0.2, 
        lwd = 0.8, arr.type = "triangle")
    shape::Arrows(arrow.pos2[1], arrow.pos2[2], arrow.pos2[3], 
        arrow.pos2[4], col = "#2f6661", arr.lwd = 1e-04, arr.length = 0.2, 
        lwd = 0.8, arr.type = "triangle")
    if (!is.null(title.name)) {
        title.pos = c(space.h/8, space.v)
        text(title.pos[1], title.pos[2], paste0(title.name, " signaling network"), 
            cex = 1)
    }
    gg <- recordPlot()
    return(gg)
}
environment(netVisual_hierarchy2) <- asNamespace('CellChat')
assignInNamespace("netVisual_hierarchy2", netVisual_hierarchy2, ns = "CellChat")

computeCommunProbPathway <- function (object = NULL, net = NULL, pairLR.use = NULL, thresh = 0.05) 
{
    if (is.null(net)) {
        net <- object@net
    }
    if (is.null(pairLR.use)) {
        pairLR.use <- object@LR$LRsig
    }
    prob <- net$prob
    prob[net$pval > thresh] <- 0
    pathways <- unique(pairLR.use$pathway_name)
    group <- factor(pairLR.use$pathway_name, levels = pathways)
    prob.pathways <- aperm(apply(prob, c(1, 2), by, group, sum), 
        c(2, 3, 1))
    pathways.sig <- pathways[apply(prob.pathways, 3, sum) != 
        0]
    prob.pathways.sig <- prob.pathways[, , pathways.sig, drop = FALSE]
    idx <- sort(apply(prob.pathways.sig, 3, sum), decreasing = TRUE, 
        index.return = TRUE)$ix
    pathways.sig <- pathways.sig[idx]
    prob.pathways.sig <- prob.pathways.sig[, , idx, drop = FALSE] ## here
    if (is.null(object)) {
        netP = list(pathways = pathways.sig, prob = prob.pathways.sig)
        return(netP)
    }
    else {
        object@netP$pathways <- pathways.sig
        object@netP$prob <- prob.pathways.sig
        return(object)
    }
}
environment(computeCommunProbPathway) <- asNamespace('CellChat')
assignInNamespace("computeCommunProbPathway", computeCommunProbPathway, ns = "CellChat")

netAnalysis_signalingRole_heatmap <- function (object, signaling = NULL, pattern = c("outgoing", "incoming", 
    "all"), slot.name = "netP", color.use = NULL, color.heatmap = "BuGn", 
    title = NULL, width = 10, height = 8, font.size = 8, font.size.title = 10, 
    cluster.rows = FALSE, cluster.cols = FALSE, show_heatmap_legend = TRUE)
{
    pattern <- match.arg(pattern)
    if (length(slot(object, slot.name)$centr) == 0) {
        stop("Please run `netAnalysis_computeCentrality` to compute the network centrality scores! ")
    }
    centr <- slot(object, slot.name)$centr
    outgoing <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
    incoming <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
    dimnames(outgoing) <- list(levels(object@idents), names(centr))
    dimnames(incoming) <- dimnames(outgoing)
    for (i in 1:length(centr)) {
        outgoing[, i] <- centr[[i]]$outdeg
        incoming[, i] <- centr[[i]]$indeg
    }
    if (pattern == "outgoing") {
        mat <- t(outgoing)
        legend.name <- "Outgoing"
    }
    else if (pattern == "incoming") {
        mat <- t(incoming)
        legend.name <- "Incoming"
    }
    else if (pattern == "all") {
        mat <- t(outgoing + incoming)
        legend.name <- "Overall"
    }
    if (is.null(title)) {
        title <- paste0(legend.name, " signaling patterns")
    }
    else {
        title <- paste0(paste0(legend.name, " signaling patterns"), 
            " - ", title)
    }
    if (!is.null(signaling)) {
        mat1 <- mat[rownames(mat) %in% signaling, , drop = FALSE]
        mat <- matrix(0, nrow = length(signaling), ncol = ncol(mat))
        idx <- match(rownames(mat1), signaling)
        mat[idx[!is.na(idx)], ] <- mat1
        dimnames(mat) <- list(signaling, colnames(mat1))
    }
    mat.ori <- mat
    mat <- sweep(mat, 1L, apply(mat, 1, max), "/", check.margin = FALSE)
    mat[mat == 0] <- NA
	mat[is.nan(mat)] <- NA # here
    if (is.null(color.use)) {
        color.use <- scPalette(length(colnames(mat)))
    }
	n_col <- if ( length(unique(as.vector(mat))) == 100 ) 101 else 100 # here, a very stupid bug in ComplexHeatmap:::ColorMapping
    color.heatmap.use = (grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, 
        name = color.heatmap))))(n_col) # here
    df <- data.frame(group = colnames(mat))
    rownames(df) <- colnames(mat)
    names(color.use) <- colnames(mat)
    col_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use), 
        which = "column", show_legend = FALSE, show_annotation_name = FALSE, 
        simple_anno_size = grid::unit(0.2, "cm"))
    ha2 = HeatmapAnnotation(Strength = anno_barplot(colSums(mat.ori), 
        border = FALSE, gp = gpar(fill = color.use, col = color.use)), 
        show_annotation_name = FALSE)
    pSum <- rowSums(mat.ori)
    pSum.original <- pSum
    pSum <- -1/log(pSum)
    pSum[is.na(pSum)] <- 0
    idx1 <- which(is.infinite(pSum) | pSum < 0)
    if (length(idx1) > 0) {
		MAX <- if ( any(pSum > 0) ) max(pSum) else max(abs(pSum))
        values.assign <- seq(MAX * 1.1, MAX * 1.5, 
            length.out = length(idx1))
        position <- sort(pSum.original[idx1], index.return = TRUE)$ix
        pSum[idx1] <- values.assign[match(1:length(idx1), position)]
    }
    ha1 = rowAnnotation(Strength = anno_barplot(pSum, border = FALSE), 
        show_annotation_name = FALSE)
    if (min(mat, na.rm = T) == max(mat, na.rm = T)) {
		legend.break <- c(0, max(mat, na.rm = T)) # here
		color.heatmap.use <- color.heatmap.use[c(1,100)] # here
		names(color.heatmap.use) <- legend.break # here
		show_heatmap_legend = FALSE # here
    }
    else {
        legend.break <- c(round(min(mat, na.rm = T), digits = 1), 
            round(max(mat, na.rm = T), digits = 1))
    }
    ht1 = Heatmap(mat, col = color.heatmap.use, na_col = "white", 
        name = "Relative strength", bottom_annotation = col_annotation, 
        top_annotation = ha2, right_annotation = ha1, cluster_rows = cluster.rows, 
        cluster_columns = cluster.rows, row_names_side = "left", 
        row_names_rot = 0, row_names_gp = gpar(fontsize = font.size), 
        column_names_gp = gpar(fontsize = font.size), width = unit(width, 
            "cm"), height = unit(height, "cm"), column_title = title, 
        column_title_gp = gpar(fontsize = font.size.title), column_names_rot = 90, 
        heatmap_legend_param = list(title_gp = gpar(fontsize = 8, 
            fontface = "plain"), title_position = "leftcenter-rot", 
            border = NA, at = legend.break, legend_height = unit(20, 
                "mm"), labels_gp = gpar(fontsize = 8), grid_width = unit(2, 
                "mm")),
		show_heatmap_legend = show_heatmap_legend) # here
    return(ht1)
}
environment(netAnalysis_signalingRole_heatmap) <- asNamespace('CellChat')
assignInNamespace("netAnalysis_signalingRole_heatmap", netAnalysis_signalingRole_heatmap, ns = "CellChat")

subsetCellChat <- function (object, cells.use = NULL, idents.use = NULL, group.by = NULL, 
    invert = FALSE, thresh = 0.05) 
{
    if (!is.null(idents.use)) {
        if (is.null(group.by)) {
            labels <- object@idents
            if (object@options$mode == "merged") {
                message("Use the joint cell labels from the merged CellChat object")
                labels <- object@idents$joint
            }
        }
        else {
            labels <- object@meta[[group.by]]
        }
        if (!is.factor(labels)) {
            labels <- factor(labels)
        }
        level.use0 <- levels(labels)
        level.use <- levels(labels)[levels(labels) %in% unique(labels)]
        if (invert) {
            level.use <- level.use[!(level.use %in% idents.use)]
        }
        else {
            level.use <- level.use[level.use %in% idents.use]
        }
        cells.use.index <- which(as.character(labels) %in% level.use)
        cells.use <- names(labels)[cells.use.index]
    }
    else if (!is.null(cells.use)) {
        labels <- object@idents
        if (object@options$mode == "merged") {
            message("Use the joint cell labels from the merged CellChat object")
            labels <- object@idents$joint
        }
        level.use0 <- levels(labels)
        level.use <- levels(labels)[levels(labels) %in% unique(as.character(labels[cells.use]))]
        cells.use.index <- which(names(labels) %in% cells.use)
    }
    else {
        stop("USER should define either `cells.use` or `idents.use`!")
    }
    cat("The subset of cell groups used for CellChat analysis are ", 
        level.use, "\n")
    if (nrow(object@data) > 0) {
        data.subset <- object@data[, cells.use.index]
    }
    else {
        data.subset <- matrix(0, nrow = 0, ncol = 0)
    }
    if (nrow(object@data.project) > 0) {
        data.project.subset <- object@data.project[, cells.use.index]
    }
    else {
        data.project.subset <- matrix(0, nrow = 0, ncol = 0)
    }
    data.signaling.subset <- object@data.signaling[, cells.use.index]
    meta.subset <- object@meta[cells.use.index, , drop = FALSE]
    if (object@options$mode == "merged") {
        idents <- object@idents[1:(length(object@idents) - 1)]
        group.existing <- level.use0[level.use0 %in% level.use]
        group.existing.index <- which(level.use0 %in% level.use)
        net.subset <- vector("list", length = length(object@net))
        netP.subset <- vector("list", length = length(object@netP))
        idents.subset <- vector("list", length = length(idents))
        names(net.subset) <- names(object@net)
        names(netP.subset) <- names(object@netP)
        names(idents.subset) <- names(object@idents[1:(length(object@idents) - 
            1)])
        images.subset <- vector("list", length = length(idents))
        names(images.subset) <- names(object@idents[1:(length(object@idents) - 
            1)])
        for (i in 1:length(idents)) {
            cat("Update slots object@images, object@net, object@netP, object@idents in dataset ", 
                names(object@idents)[i], "\n")
            images <- object@images[[i]]
            for (images.j in names(images)) {
                values <- images[[images.j]]
                if (images.j %in% c("coordinates")) {
                  values.new <- values[cells.use.index, ]
                  images[[images.j]] <- values.new
                }
                if (images.j %in% c("distance")) {
                  values.new <- values[group.existing.index, 
                    group.existing.index, drop = FALSE]
                  images[[images.j]] <- values.new
                }
            }
            images.subset[[i]] <- images
            net <- object@net[[i]]
            for (net.j in names(net)) {
                values <- net[[net.j]]
                if (net.j %in% c("prob", "pval")) {
                  values.new <- values[group.existing.index, 
                    group.existing.index, ]
                  net[[net.j]] <- values.new
                }
                if (net.j %in% c("count", "sum", "weight")) {
                  values.new <- values[group.existing.index, 
                    group.existing.index]
                  net[[net.j]] <- values.new
                }
            }
            net.subset[[i]] <- net
            netP = computeCommunProbPathway(net = net.subset[[i]], 
                pairLR.use = object@LR[[i]]$LRsig, thresh = thresh)
#            netP$centr = netAnalysis_computeCentrality(net = net.subset[[i]]$prob) ## ???
            netP$centr = netAnalysis_computeCentrality(net = netP$prob) ## here
            netP.subset[[i]] <- netP
            idents.subset[[i]] <- idents[[i]][names(idents[[i]]) %in% 
                cells.use]
            idents.subset[[i]] <- factor(idents.subset[[i]], 
                levels = levels(idents[[i]])[levels(idents[[i]]) %in% 
                  level.use])
        }
        idents.subset$joint <- factor(object@idents$joint[cells.use.index], 
            levels = level.use)
    }
    else {
        cat("Update slots object@images, object@net, object@netP in a single dataset...", 
            "\n")
        group.existing <- level.use0[level.use0 %in% level.use]
        group.existing.index <- which(level.use0 %in% level.use)
        images <- object@images
        for (images.j in names(images)) {
            values <- images[[images.j]]
            if (images.j %in% c("coordinates")) {
                values.new <- values[cells.use.index, ]
                images[[images.j]] <- values.new
            }
            if (images.j %in% c("distance")) {
                values.new <- values[group.existing.index, group.existing.index, 
                  drop = FALSE]
                images[[images.j]] <- values.new
            }
        }
        images.subset <- images
        net <- object@net
        for (net.j in names(net)) {
            values <- net[[net.j]]
            if (net.j %in% c("prob", "pval")) {
                values.new <- values[group.existing.index, group.existing.index, ,
                  drop = FALSE]
                net[[net.j]] <- values.new
            }
            if (net.j %in% c("count", "sum", "weight")) {
                values.new <- values[group.existing.index, group.existing.index, 
                  drop = FALSE]
                net[[net.j]] <- values.new
            }
        }
        net.subset <- net
        netP = computeCommunProbPathway(net = net.subset, pairLR.use = object@LR$LRsig, 
            thresh = thresh)
#        netP$centr = netAnalysis_computeCentrality(net = net.subset$prob) ## ???
        netP$centr = netAnalysis_computeCentrality(net = netP$prob) ## here
        netP.subset <- netP
        idents.subset <- object@idents[cells.use.index]
        idents.subset <- factor(idents.subset, levels = level.use)
    }
    object.subset <- methods::new(Class = "CellChat", data = data.subset, 
        data.signaling = data.signaling.subset, data.project = data.project.subset, 
        images = images.subset, net = net.subset, netP = netP.subset, 
        meta = meta.subset, idents = idents.subset, var.features = object@var.features, 
        LR = object@LR, DB = object@DB, options = object@options)
    return(object.subset)
}
environment(subsetCellChat) <- asNamespace('CellChat')
assignInNamespace("subsetCellChat", subsetCellChat, ns = "CellChat")


