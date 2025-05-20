Software :
	R-4.1.3
	CellChat-1.6.1

Script : 
	Step1_Seurat2CellChat.sh : 
		1. Read DataBase
		2. Convert Seurat object to CellChat object

	Step2_Run_CellChat.sh :
		1. identifyOverExpressedInteractions
		2. computeCommunProb
		3. aggregateNet (thresh = 0.05)
		4. computeCommunProbPathway
		5. draw figure (network, dotplot, contribution)
