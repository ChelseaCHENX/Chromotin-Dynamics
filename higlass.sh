# https://github.com/higlass/higlass/wiki/Setup#running-locally
pip install higlass-manage
pip install higlass-python

# does not allow links between mac and box drive => cp mcools to Desktop/mcools
higlass-manage ingest SRR6493702.GRCh38.mapq_30.50000.mcool --project-name mcf7
mcools  for f in *; do higlass-manage ingest $f --project-name mcf7; done