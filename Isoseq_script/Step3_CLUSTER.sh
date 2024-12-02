## cluster
isoseq3 cluster --log-level INFO 2.FLNC/flnc.consensusreadset.xml 3.Cluster/unpolished.transcriptset.xml --use-qvs -j 4
isoseq3 summarize --log-level INFO 3.Cluster/unpolished.transcriptset.xml primers.fasta 3.Cluster/summary.csv

## polish
samtools-1.9 view -h 3.Cluster/unpolished.hq.bam | perl -ne 'if ( ! /^@/ ) { my @data = split /\t/; $data[10] = "~" x length $data[9]; $_ = join("\t", @data); } print;' | samtools-1.9 view -Sb - > 4.Isoform/hq_transcripts.bam
samtools-1.9 view -h 3.Cluster/unpolished.lq.bam | perl -ne 'if ( ! /^@/ ) { my @data = split /\t/; $data[10] = "~" x length $data[9]; $_ = join("\t", @data); } print;' | samtools-1.9 view -Sb - > 4.Isoform/lq_transcripts.bam
pbindex 4.Isoform/hq_transcripts.bam
pbindex 4.Isoform/lq_transcripts.bam
dataset create --type TranscriptSet 4.Isoform/hq.transcriptset.xml 4.Isoform/hq_transcripts.bam --force
dataset create --type TranscriptSet 4.Isoform/lq.transcriptset.xml 4.Isoform/lq_transcripts.bam --force
samtools-1.9 view 4.Isoform/hq_transcripts.bam | awk '{print ">"$1"\n"$10}' > 4.Isoform/hq_transcripts.fasta
samtools-1.9 view 4.Isoform/lq_transcripts.bam | awk '{print ">"$1"\n"$10}' > 4.Isoform/lq_transcripts.fasta
samtools-1.9 view 4.Isoform/hq_transcripts.bam | awk '{print "@"$1"\n"$10"\n+\n"$11}' > 4.Isoform/hq_transcripts.fastq
samtools-1.9 view 4.Isoform/lq_transcripts.bam | awk '{print "@"$1"\n"$10"\n+\n"$11}' > 4.Isoform/lq_transcripts.fastq
python3 -m pbreports.report.isoseq3 --log-level INFO 4.Isoform/hq.transcriptset.xml 4.Isoform/lq.transcriptset.xml 4.Isoform/isoseq3.report.json
