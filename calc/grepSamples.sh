# TCC
awk 'FS="\t", OFS="\t" { gsub("ftp.sra.ebi.ac.uk", "era-fasp@fasp.sra.ebi.ac.uk:"); print }' PRJNA428600 | cut -f 11| awk -F ";" 'OFS="\n" {print $1, $2}' | \
awk NF | awk 'NR > 1, OFS="\n" {print "ascp -QT -l 300m -P33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh" " " $1 " ."}' > download.txt

# atac-seq https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144580
bash
for x in {10992898..10992915}; do prefetch SRR$x &;done
for x in {10992902..10992915}; do vdb-validate SRR$x ;done 2&>1 >> checksum.txt 
for x in {10992902..10992915}; do fastq-dump --split-files SRR$x/SRR$x.sra &;done 
for x in {10992898..10992901}; do fastq-dump --split-files SRR$x/SRR$x.sra &;done 

for x in {10992898..10992915}
do fastq-dump --split-files SRR$x --gzip >> SRR$x.download.log 2>&1 &
done
