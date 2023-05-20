#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, sys, getopt
import pandas as pd
from Bio import Entrez 

#================================================
# FUNCTION get_taxonomy(Organism)
# INPUT : Organism = organism name
# OUTPUT : organism limeage (dictionary)
# get complete taxonomy from an organism name using NCBI EUtilities
# code is customisable to get different levels of taxonomy
#================================================
def get_taxonomy(Organism):
    # log in to have a limit of 10 requests / second instead of 3
    # login e-mail should be modified
    Entrez.email = "mathe.meije@gmail.com"
    # search for the Organism name in the NCBI taxonomy database
    handle = Entrez.esearch(db = "taxonomy", term =  Organism, retmax = 100)
    print(Organism)
    records = Entrez.read(handle)
    # retrieve the organism NCBI Id
    identifiers = records['IdList']
    taxonomy = {}
    if identifiers :
        # use the NCBI Id to fetch information from NCBI (xml format)
        handle = Entrez.efetch(db = "taxonomy", id = identifiers, retmax = 100)
        records = Entrez.read(handle)
        # read xml file and extract lineage, stored in a dict (taxonomy)
        lineageEx = records[0]['LineageEx']
        for lineage in lineageEx:
            if lineage['Rank'] == 'superkingdom':
                taxonomy['Superkingdom'] = lineage['ScientificName']
            if lineage['Rank'] == 'kingdom':
                taxonomy['Kingdom'] = lineage['ScientificName']
            if lineage['Rank'] == 'phylum':
                taxonomy['Phylum'] = lineage['ScientificName']
            if lineage['Rank'] == 'class':
                taxonomy['Class'] = lineage['ScientificName']
            if lineage['Rank'] == 'order':
                taxonomy['Order'] = lineage['ScientificName']
            if lineage['Rank'] == 'family':
                taxonomy['Family'] = lineage['ScientificName']
            if lineage['Rank'] == 'genus':
                taxonomy['Genus'] = lineage['ScientificName']
    else :
        taxonomy['Superkingdom'] = 'None'
    return(taxonomy)

#================================================
# FUNCTION filter_TABLE_res(TABLE_res, coverage, identity)
# INPUT : TABLE_res = dataframe containing results from BLAST 
# PARAMS : coverage = coverage threshold for filtering
#          identity = identity threshold for filtering
# OUTPUT : filtered dataframe (coverage, threshold, order and family)
# filter the table based on coverage, identity and taxonomy 
#================================================
def filter_TABLE_res(TABLE_res, coverage, identity):
    # compute coverage in new dataframe column 
    TABLE_res['Coverage'] = (TABLE_res['Length'] / TABLE_res['queryLength'])*100
    # apply coverage and identity thresholds
    TABLE_filtered = TABLE_res.query('Coverage >=' + str(coverage) + '& Identities >=' + str(identity))
    # get complete taxonomy for all organisms in the filtered table
    for ind in TABLE_filtered.index:
        taxonomy = get_taxonomy(TABLE_filtered['Organism'][ind])
        tax_keys = ['Superkingdom', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus']
        for k in tax_keys :
            if k in taxonomy :
                TABLE_filtered.at[ind, k] = taxonomy[k]
            else :
                TABLE_filtered.at[ind, k] = float('nan')
    # remove rows if the family is unknown
    TABLE_filtered.dropna(subset = ['Family'], inplace = True)
    # sort rows : order and family alphabetically, ascending coverage
    TABLE_filtered_sorted = TABLE_filtered.sort_values(['Order', 'Family', 'Coverage'], ascending=(True, True, True))
    # keep one representtative / family, the one with the highest coverage
    TABLE_filtered_uniq_sorted_fam = TABLE_filtered_sorted.drop_duplicates(subset = ['Family'], keep = "last")
    return(TABLE_filtered_uniq_sorted_fam)

#================================================
# FUNCTION read_table(FILIN, queryLength = 0)
# INPUT : tsv file (BLAST result)
# PARAMS : BLAST query sequence length
# OUTPUT : pandas dataframe
#================================================
def read_table(FILIN, queryLength):
    # load input file (tsv format) in a dataframe
    TABLE_res = pd.read_csv(FILIN, sep='\t')
    # rename dataframe columns
    TABLE_res.columns = ['Hit', 'DB', 'Accession', 'Description', 'Organism', 'Length', 'Score', 'Identities', 'Positives', 'E']
    # add the query sequence length in the table
    TABLE_res['queryLength'] = [int(queryLength)] * 1000
    pd.to_numeric(TABLE_res['queryLength'])
    # add empty columns for taxonomy infos
    TABLE_res["Superkingdom"] = ""
    TABLE_res["Kingdom"] = ""
    TABLE_res["Phylum"] = ""
    TABLE_res["Class"] = ""
    TABLE_res["Order"] = ""
    TABLE_res["Family"] = ""
    TABLE_res["Genus"] = ""
    return(TABLE_res)

#================================================
# 
#================================================
def main(argv):
    opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
    for opt, arg in opts:
        if opt == '-h':
            print ('test.py -i <inputfile> -o <outputfile>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            FILIN = arg
        elif opt in ("-o", "--ofile"):
            FILOUT = arg
    if not os.path.isfile(FILIN):
        print(FILIN + " does not exists.")
        exit()
    # please provide query sequence length
    # 607 input length
    TABLE = read_table(FILIN, queryLength = 607)
    # please provide coverage and identity thresholds
    TABLE_filtered = filter_TABLE_res(TABLE, coverage = 90, identity = 90)
    print(TABLE_filtered)
    TABLE_filtered.to_csv(FILOUT, sep="\t")

if __name__ == "__main__":
    main(sys.argv[1:])