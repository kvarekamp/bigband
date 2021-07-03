# Bigband - BigBird trained on DNA sequences

## genbank folder
This folder contains the genbank database. It can be copied over ftp from https://ftp.ncbi.nlm.nih.gov/genbank/

There is a lot of data there (about 1/2 a TB) in the form of hundreds of gzips. The git repo contains one (gbenv18.seq.gz) to get started - the rest should be downloaded. The gzips contain flat ascii files which are human readable - unzip one to have a look.

There is also a README.genbank in there that explains the file format in more detail.

## read_genbank.py 
This script iterates through the genbank folder, reads all the gzips in there and parses them into dictionaries.

Just run it to see an example of the output in dump.json


