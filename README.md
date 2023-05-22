## NCBI

### Content
**BLAST_filter.py** : filter a BLAST result table using coverage, identity and taxonomy.  
**download_NCBI_genomes.py** : download GenBank genomes from NCBI ftp site using a list of GenBank IDs.  
**get_NCBI_metadata.py** : retrieve (customizable) metadata from a set of BioSample IDs using NCBI EUtilities and BioPython. 

### Usage
- BLAST_filter.py -i <inputfile> -o <outdir>
```bash
./BLAST_filter.py -i BLAST_result_table.tsv -o output_file.tsv
```
- download_NCBI_genomes.py <inputfile> <outputfile>
```bash
./download_NCBI_genomes.py genomesID.txt genomesDirectory/
# unzip genome files
cd genomesDirectory/
gunzip *
```
Input file (txt format) : list of GenBank IDs
```bash
GCA_002686955.1  
GCA_002898075.1  
GCA_013042495.1  
GCA_001795235.1  
GCA_003229755.1
GCA_003349885.1
GCA_003230085.1
GCA_011371465.1
GCA_003696725.1  
```
