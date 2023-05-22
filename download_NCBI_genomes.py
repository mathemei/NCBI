#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os, requests, re

def main(file, outDir):
	FILE_NAME = file
	PATH_OUT = outDir
	if not os.path.isfile(FILE_NAME):
		print(FILE_NAME + "does not exists.")
		exit()
	# get GenBank ids list
	with open(FILE_NAME) as f:
		GenBankIDs = f.read().splitlines()
	# remove N/A
	GenBankIDs = [i for i in GenBankIDs if i != "N/A"]
	# only if one line = ID + taxonomy
	GenBankIDs = [i.split(",")[0] for i in GenBankIDs]

	for ID in GenBankIDs:
		url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/" + ID[0:3] + "/" + ID[4:7] + "/" + ID[7:10] + "/" + ID[10:13] + "/" 
		r = requests.get(url, allow_redirects=True)
		TMP_FILE = "tmp" 
		open(TMP_FILE, 'wb').write(r.content)
		textfile = open(TMP_FILE, 'r')
		matches = []
		reg = re.compile("\"" + ID + ".+\"")
		for line in textfile:
			matches += reg.findall(line)
		textfile.close()
		os.remove(TMP_FILE)
		url = url + matches[0][1:-2] + "/" + matches[0][1:-2] + "_genomic.fna.gz"
		GENOME_FILE = PATH_OUT + ID + ".fna.gz"
		r = requests.get(url, allow_redirects=True)
		open(GENOME_FILE, 'wb').write(r.content)

if __name__ == "__main__":
	main(sys.argv[1], sys.argv[2])
