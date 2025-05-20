lima --isoseq --peek-guess -j 1 --alarms 2.FLNC/alarms.json 1.CCS/ccs.bam primers.fasta 2.FLNC/fl.datastore.json
python3 -c 'import sys; import os.path ; from pbcommand.models import DataStore ; ds = DataStore.load_from_json(os.path.realpath(sys.argv[1])) ; print("\n".join([f.path for f in ds.files.values()]))' 2.FLNC/fl.datastore.json > 2.FLNC/datasets.fofn
dataset --strict create --unique-collections --no-sub-datasets 2.FLNC/merged.consensusreadset.xml 2.FLNC/datasets.fofn --force
isoseq3 refine --log-level INFO --require-polya -j 4 2.FLNC/merged.consensusreadset.xml /Bio/Bin/pipeline/IsoSeq/isoseq_shell/dev/lib/primers.fasta 2.FLNC/flnc.consensusreadset.xml 2.FLNC/flnc.filter_summary.json 2.FLNC/flnc.report.csv
perl -MJSON -e 'my $data = JSON->new->decode(join "", <>); my %out = map { $_->{id}, $_->{value} } @{$data->{attributes}}; print to_json(\%out, {pretty=>1});' 2.FLNC/flnc.filter_summary.json > 2.FLNC/flnc.filter_summary.oldstyle.json
python3 -m pbreports.report.barcode_isoseq3 --log-level INFO 2.FLNC/fl.datastore.json 1.CCS/gathered.consensusreadset.xml primers.fasta 2.FLNC/barcode_isoseq3.report.json 2.FLNC/barcode_isoseq3.report.csv 2.FLNC/flnc.filter_summary.oldstyle.json 2.FLNC/flnc.consensusreadset.xml
rm -f 2.FLNC/flnc.filter_summary.oldstyle.json
samtools-1.9 view 2.FLNC/flnc.bam | awk '{print ">"$1"\n"$10}' > 2.FLNC/flnc.fasta

