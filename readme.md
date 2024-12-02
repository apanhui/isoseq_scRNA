### software :
- smrtlink_9.0.0.92188
- samtools-1.9
- cd-hit-v4.6.7
- diamond-0.9.25
- R-4.1.3
- DoubletFinder-2.0.3
- Seurat-4.1.0

### Isoseq script :
- Isoseq_script/Step1_CCS.sh     : Consensus generation
- Isoseq_script/Step2_FLNC.sh    : Primer Removal
- Isoseq_script/Step3_CLUSTER.sh : Isoseq3 cluster & Polish
- Isoseq_script/Step4_CDHIT.sh   : Cluster hq_transcripts
- Isoseq_script/Step5_Annot.sh   : Gene Annotation

### Seurat script :
- Seurat_script/Step1_DoubletFinder.sh : Remove Doublet
- Seurat_script/Step2_Seurat.sh        : Run Seurat & draw figure

