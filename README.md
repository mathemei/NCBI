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
- get_NCBI_metadata.py <inputfile> 
```bash
./BLAST_filter.py genomesFile.tsv
```
Input file (tsv format) : table from NCBI Taxonomy search page  
```bash
Assembly Accession	Assembly Name	Organism Name	Organism Infraspecific Names Breed	Organism Infraspecific Names Strain	Organism Infraspecific Names Cultivar	Organism Infraspecific Names Ecotype	Organism Infraspecific Names Isolate	Organism Infraspecific Names Sex	Annotation Name	Assembly Stats Total Sequence Length	Assembly Level	Assembly Submission Date	WGS project accession	Assembly BioProject Accession	Assembly BioSample Accession
GCF_000172135.1	ASM17213v1	Bifidobacterium dentium ATCC 27678		ATCC 27678					NCBI Prokaryotic Genome Annotation Pipeline (PGAP)	2642081	Contig	2008-03-13	ABIX02	PRJNA20555	SAMN00000021
GCF_000196555.1	ASM19655v1	Bifidobacterium longum subsp. longum JCM 1217		JCM 1217					NCBI Prokaryotic Genome Annotation Pipeline (PGAP)	2385164	Complete Genome	2011-01-28		PRJDA32047	SAMD00060951
GCF_000224965.2	ASM22496v2	Bifidobacterium animalis subsp. lactis BLC1		BLC1					NCBI Prokaryotic Genome Annotation Pipeline (PGAP)	1938583	Complete Genome	2013-02-05		PRJNA71815	SAMN02604239
GCF_000273525.1	PB_Bifi_bifi_NCIMB_41171_V1	Bifidobacterium bifidum NCIMB 41171		NCIMB 41171					NCBI Prokaryotic Genome Annotation Pipeline (PGAP)	2216395	Contig	2012-06-15	AKCA01	PRJNA30055	SAMN02463676
GCF_000304215.1	ASM30421v1	Bifidobacterium asteroides PRL2011		PRL2011					NCBI Prokaryotic Genome Annotation Pipeline (PGAP)	2167304	Complete Genome	2012-10-09		PRJNA83029	SAMN02604114
```

### Bibliography
