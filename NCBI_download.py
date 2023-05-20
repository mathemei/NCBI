#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys, os, requests, re
import pandas as pd

# input : tsv format, downloaded from NCBI genomes search results (custom columns)
# columns : Assembly Accession - Assembly Name - Organism Name - Organism Infraspecific Names Breed    - Organism Infraspecific Names Strain - 
# Organism Infraspecific Names Cultivar    - Organism Infraspecific Names Ecotype - Organism Infraspecific Names Isolate - Organism Infraspecific Names Sex    
# Annotation Name - Assembly Stats Total Sequence Length - Assembly Level - Assembly Submission Date - WGS project accession - Assembly BioProject Accession - Assembly BioSample Accession

#================================================
# INPUT : tsv file from NCBI Genomes search results
# downloads all genomes using their assembly accession 
#================================================
def main(genomeFile) :
    # test if input file exists
    if not os.path.isfile(genomeFile):
        print(genomeFile + " does not exists.")
        exit()
    # output directory, custom
    PATH_OUT = "./faecalibacterium/"
    if not os.path.exists(PATH_OUT):
        print("Creating output directory : " + PATH_OUT)
        os.makedirs(PATH_OUT)
    # read input file in pandas dataframe
    genomes = pd.read_csv(genomeFile, sep='\t')
    # retrieve genomes accession ID
    genomeIds = genomes["Assembly Accession"]
    for ID in genomeIds :
        # create URL from NCBI ftp
        url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/" + ID[0:3] + "/" + ID[4:7] + "/" + ID[7:10] + "/" + ID[10:13] + "/" 
        r = requests.get(url, allow_redirects=True)
        TMP_FILE = "tmp" 
        # temp file with the content of ftp page corresponding to the genome ID
        open(TMP_FILE, 'wb').write(r.content)
        # retrieve the assembly version from temp file
        textfile = open(TMP_FILE, 'r')
        matches = []
        reg = re.compile("\"" + ID + ".+\"")
        for line in textfile:
            matches += reg.findall(line)
        textfile.close()
        os.remove(TMP_FILE)
        # get the URL to the genome file
        url = url + matches[0][1:-2] + "/" + matches[0][1:-2] + "_genomic.fna.gz"
        # download genome
        GENOME_FILE = PATH_OUT + ID + ".fna.gz"
        r = requests.get(url, allow_redirects=True)
        open(GENOME_FILE, 'wb').write(r.content)

if __name__ == "__main__":
    main(sys.argv[1])
