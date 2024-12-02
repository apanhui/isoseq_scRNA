## ccs chunk
mkdir -p 1.CCS/subccs1
perl -e 'for my $i ( 1 .. 100 ) { my $path = "1.CCS/subccs1/chunk$i"; print "mkdir -p $path; ccs-alt --log-level INFO m64269e_230627_085133.subreads.bam $path/out.consensusreadset.xml --chunk $i/100 --min-length 50 --max-length 15000 --min-passes 1 --min-snr 2.5 --min-rq 0.8 -j 4; samtools-1.9 index $path/out.bam\n"; }' > 1.CCS/subccs1/run.sh
sh 1.CCS/subccs1/run.sh >> 1.CCS/subccs1/run.sh.log 2>&1

if [[ ! `ls -1 1.CCS/subccs1/chunk*/out.consensusreadset.xml | wc -l` -eq 100 ]]; then exit 1; fi
rm -rf 1.CCS/subccs2
mkdir -p 1.CCS/subccs2
perl -e 'for my $i ( 1 .. 100 ) { my $path = "1.CCS/subccs2/chunk$i"; print "mkdir -p $path; ccs-alt --log-level INFO FISO23H001595-1A/m64268e_230707_092127.subreads.bam $path/out.consensusreadset.xml --chunk $i/100 --min-length 50 --max-length 15000 --min-passes 1 --min-snr 2.5 --min-rq 0.8 -j 4; samtools-1.9 index $path/out.bam\n"; }' > 1.CCS/subccs2/run.sh
sh 1.CCS/subccs2/run.sh >> 1.CCS/subccs2/run.sh.log 2>&1

if [[ ! `ls -1 1.CCS/subccs2/chunk*/out.consensusreadset.xml | wc -l` -eq 100 ]]; then exit 1; fi
rm -rf 1.CCS/subccs3
mkdir -p 1.CCS/subccs3
perl -e 'for my $i ( 1 .. 100 ) { my $path = "1.CCS/subccs3/chunk$i"; print "mkdir -p $path; ccs-alt --log-level INFO m64268e_230707_092127.subreads.bam $path/out.consensusreadset.xml --chunk $i/100 --min-length 50 --max-length 15000 --min-passes 1 --min-snr 2.5 --min-rq 0.8 -j 4; samtools-1.9 index $path/out.bam\n"; }' > 1.CCS/subccs3/run.sh
sh 1.CCS/subccs3/run.sh >> 1.CCS/subccs3/run.sh.log 2>&1

if [[ ! `ls -1 1.CCS/subccs3/chunk*/out.consensusreadset.xml | wc -l` -eq 100 ]]; then exit 1; fi

## ccs post
dataset --log-level DEBUG --strict create --name 'CCS Data' --unique-collections 1.CCS/gathered.consensusreadset.xml 1.CCS/*/chunk*/out.consensusreadset.xml --force
dataset --log-level INFO consolidate 1.CCS/gathered.consensusreadset.xml 1.CCS/ccs.bam 1.CCS/ccs.xml
python3 -m pbreports.report.ccs --log-level INFO --output-dir 1.CCS 1.CCS/gathered.consensusreadset.xml 1.CCS/ccs.report.json 1.CCS/ccs.report.csv
python3 -m pbcoretools.tasks.bam2fastq_archive --log-level DEBUG 1.CCS/gathered.consensusreadset.xml 1.CCS/ccs.fastq.zip
/Bio/bin/samtools-1.9 view 1.CCS/ccs.bam | awk '{print ">"$1"\n"$10}' > 1.CCS/ccs.fasta

