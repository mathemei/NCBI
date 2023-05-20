#!/usr/bin/env python3
# -*- coding: utf-8 -*

import sys, os, time, csv, math
import pandas as pd
from Bio import Entrez 
import xml.etree.ElementTree as ET
import numpy as np

#================================================
# FUNCTION get_metadata(BioSampleID, query)
# INPUT : BioSampleID = NCBI BioSample ID
# PARAMS : query = type of metadata to retrieve
# OUTPUT : metadata requiered from a given sample (if it exists)
#================================================
def get_metadata(BioSampleID, query):
    # log in to have a limit of 10 requests / second instead of 3
    # login e-mail should be modified
    Entrez.email = "mathe.meije@gmail.com"
    # get xml content of the Biosample entry
    handle = Entrez.efetch(db = "biosample", id = BioSampleID, retmax = 100, rettype = "xml", retmode = "xml")
    # parse xml file with ElementTree
    records = ET.parse(handle)
    BioSampleSet = records.getroot()
    for BioSample in BioSampleSet :
        for B in BioSample :
            if B.tag == "Attributes" :
                for Attribute in B :
                    # extract field corresponding to the wanted metadata
                    if 'harmonized_name' in Attribute.attrib and Attribute.attrib['harmonized_name'] == query:
                        return(Attribute.text)
#================================================
#
#================================================
def main(genomeFile):
    # test if input file exists
    if not os.path.isfile(genomeFile):
        print(genomeFile + " does not exists.")
        exit()
    # read file content in pandas dataframe
    genomes = pd.read_csv(genomeFile, sep='\t')
    # extract BioSample NCBI IDs
    BioSampleIds = genomes["Assembly BioSample Accession"]
    # handling of the EUtilities random HTTP Error 400
    # verify if a previous metadata file exists
    if not os.path.isfile("metadata_" + genomeFile):
        genomes['Isolation Source'] = "NA"
    else :
        metadata = pd.read_csv("metadata_" + genomeFile, sep='\t')
        genomes['Isolation Source'] = metadata['Isolation Source']
    genomes = genomes.fillna('NA')
    # retrieve metadata
    for ID in BioSampleIds :
        print(ID)
        print(genomes.at[genomes.index[genomes['Assembly BioSample Accession'] == ID][0], 'Isolation Source'])
        # if metadata was not previously retrieved, fetch it
        if genomes.at[genomes.index[genomes['Assembly BioSample Accession'] == ID][0], 'Isolation Source'] == "NA" :
            IS = get_metadata(ID, "isolation_source")
            print(IS)
            # case where metadata doesn't exist for the sample
            if IS == None :
                genomes.at[genomes.index[genomes['Assembly BioSample Accession'] == ID][0], 'Isolation Source'] = "missing"
            else:
                genomes.at[genomes.index[genomes['Assembly BioSample Accession'] == ID][0], 'Isolation Source'] = IS
            # write metadata file for each sample
            genomes.to_csv("metadata_" + genomeFile, sep = '\t')
        # helps with the HTTP Error but I don't really know why (reduce nb of requests / seconds maybe ?)
        time.sleep(1)

if __name__ == "__main__":
    main(sys.argv[1])