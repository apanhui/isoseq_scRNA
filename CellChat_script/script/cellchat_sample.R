

if ( sys.nframe() == 0 ) {
		library(CellChat)

		bin <- dirname(normalizePath(sub('--file=', '',  grep('--file=', commandArgs(), value = T))))
		source(paste0(bin, "/cellchat_lib.R"))
		source(paste0(bin, "/modify_cellchat.R"))

		args <- commandArgs(T)
		infile <- args[1]
		outdir <- args[2]

		cellchat <- readRDS(infile)

		dir.create(outdir, F, T)
		setwd(outdir)
}

#cellchat <- TenX2CellChat("/state/partition1/WORK2/Bio/Project/xsy_filter/singlecell/rawdata/RNA/10k_PBMC_3p_25-30yr_female")
#future::plan("multicore", workers = 4)
#cellchat <- identifyOverExpressedGenes(cellchat, thresh.p = 2, only.pos = F, thresh.fc = -1, thresh.pc = -1)
#cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat, features = rownames(cellchat@data.signaling))

#cellchat <- projectData(cellchat, PPI.human)

if ( nrow(cellchat@LR$LRsig) == 0 ) {
		stop("No LRsig found",
			"\nnrow(DB) = ", nrow(cellchat@DB$interaction),
			"\nnorw(data.signalin) = ", nrow(cellchat@data.signaling))
}

##### I. Inference of cell-cell communication network
### 1. Compute the communication probability 
cellchat@idents <- droplevels(cellchat@idents)
writeLines(levels(cellchat@idents), 'idents.list')

population.size <- if ( do.call(`/`, as.list(rev(range(table(cellchat@idents))))) > 20 ) TRUE else FALSE
cellchat <- computeCommunProb(cellchat, population.size = population.size)
cellchat <- SelectCommunication(cellchat,
		source = cellchat@options$select_commun$source,
		target = cellchat@options$select_commun$target,
		method = cellchat@options$select_commun$method)
cellchat <- aggregateNet(cellchat, thresh = 0.05)

StatCommunProb(cellchat, thresh = 0.05)
PlotNet(cellchat)
PlotNet(cellchat, split.group = T)
PlotLRDotPlot(cellchat, thresh = 0.05)

data <- GetOnlineData(cellchat)
WriteTable(data, file = 'CellChat.LR.all.xls')
data <- subset(data, Pvalue <= 0.05)
write.table(data, file = 'CellChat.LR.filtered.xls', row.names = F, sep = "\t", quote = F)


### 2. Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat, thresh = 0.05)
data <- StatCommunProb(cellchat, outpref = "Pathway.CommunProb", slot.name = "netP", thresh = 0.05)

pathways.show <- (subset(data, prob > 0) %>% group_by(pathway_name) %>% summarise(pval = sum(pval), prob = sum(prob)) %>% arrange(pval, -prob) %>% head(5))[['pathway_name']]
writeLines(pathways.show, 'pathways.show.list')

for ( i in pathways.show ) {
		pdf(paste0("Pathway.", i, ".hierarchy.pdf"), 6, 6)
		netVisual_aggregate(cellchat, signaling = i, layout = 'hierarchy', vertex.receiver = 1:round(nlevels(cellchat@idents) / 2))
		dev.off()

		pdf(paste0("Pathway.", i, ".networks.pdf"), 6, 6)
		netVisual_aggregate(cellchat, signaling = i, layout = "circle")
		dev.off()

		w <- max(sapply(unique(cellchat@idents), function(i) grid::convertWidth(grid::stringWidth(i), 'in')))
		pdf(paste0("Pathway.", i, ".circle.pdf"), 6 + w*2, 6 + w*2)
		netVisual_aggregate(cellchat, signaling = i, layout = "chord")
		dev.off()

		pdf(paste0("Pathway.", i, ".heatmap.pdf"), 6, 6)
		ht <- netVisual_heatmap(cellchat, signaling = i, color.heatmap = c("white", "darkred"))
#		draw(ht)
		print(ht)
		dev.off()

}

df <- NULL
for ( i in cellchat@netP$pathways ) {
		p <- netAnalysis_contribution(cellchat, signaling = i, return.data = T)
		tmp <- p$LR.contribution
		tmp$pathway_name <- i
		df <- rbind(df, tmp) 
}
df <- df %>% select(pathway_name, Pair = name, contribution)
write.table(df, file = "Pathway.contribution.xls", row.names = F, sep = "\t", quote = F)

for ( i in pathways.show ) {
		p <- netAnalysis_contribution(cellchat, signaling = i)
		ggsave(p, file = paste0("Pathway.", i, ".contribution.pdf"), width = 6, height = 6)

		p <- netVisual_bubble(cellchat, signaling = i, remove.isolate = T)
		wh <- GetWH(p, default.height = 4)
		ggsave(p, file = paste0("Pathway.", i, ".dotplot.pdf"), width = wh[1] + 1, height = wh[2])
}




############### below result won't be shown in report ################

##### II. Systems analysis of cell-cell communication network
### 1. Identify signaling roles of cell groups as well as the major contributing signaling
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 

for ( i in pathways.show ) {
		pdf(paste0("SignalingRole.", i, ".network.pdf"), 6, 6)
		ht <- netAnalysis_signalingRole_network(cellchat, signaling = i, width = 8, height = 2.5, font.size = 10)
		print(ht)
		dev.off()
}

gg1 <- netAnalysis_signalingRole_scatter(cellchat, do.label = T)
ggsave(gg1, file = "SignalingRole.scatter.pdf", width = 6, height = 6)
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = pathways.show, do.label = T)
ggsave(gg2, file = "SignalingRole.scatter.top5_pathway.pdf", width = 6, height = 6)

for ( i in pathways.show) {
		p <- netAnalysis_signalingRole_scatter(cellchat, signaling = i, do.label = F)
		p <- p + guides(colour = guide_legend(title = 'Cluster', order = 2))
		ggsave(p, file = paste0("Pathway.", i, ".scatter.pdf"), width = 6, height = 6)
}

ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
pdf("SignalingRole.heatmap.pdf", 12, 6)
ht1 + ht2
dev.off()

### 2. Identify global communication patterns to explore how multiple cell types and signaling pathways coordinate together
if ( 0 ) { ## unable automatize yet
library(NMF)
library(ggalluvial)

p <- selectK(cellchat, pattern = "outgoing")
ggsave(p, file = "NMF.outgoing.k.select.pdf",width = 8, height = 6)

nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)

p <- netAnalysis_river(cellchat, pattern = "outgoing")
ggsave(p, file = "SignalingPattern.Outgoing.river.pdf", width = 7, height = 6)

p <- netAnalysis_dot(cellchat, pattern = "outgoing")
ggsave(p, file = "SignalingPattern.Outgoing.dot.pdf", width = 7, height = 6)


p <- selectK(cellchat, pattern = "incoming")
ggsave(p, file = "NMF.incoming.k.select.pdf",width = 8, height = 6)

nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)

p <- netAnalysis_river(cellchat, pattern = "incoming")
ggsave(p, file = "SignalingPattern.Incoming.river.pdf", width = 7, height = 6)

p <- netAnalysis_dot(cellchat, pattern = "incoming")
ggsave(p, file = "SignalingPattern.Incoming.dot.pdf", width = 7, height = 6)

}

### 3. Manifold and classification learning analysis of signaling networks
if ( dim(cellchat@netP$prob)[3] > 2 ) {
cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional", umap.method = 'uwot')
cellchat <- netClustering(cellchat, type = "functional")
p <- netVisual_embedding(cellchat, type = "functional", label.size = 3.5)
ggsave(p, file = "NetSimilarity.functional.pdf", width = 6, height = 6)


cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural", umap.method = 'uwot')
cellchat <- netClustering(cellchat, type = "structural")
p <- netVisual_embedding(cellchat, type = "structural", label.size = 3.5)
ggsave(p, file = "NetSimilarity.structural.pdf", width = 6, height = 6)
}

##### END #####
saveRDS(cellchat, file = 'cellchat.Rds')

message('_complete_')



