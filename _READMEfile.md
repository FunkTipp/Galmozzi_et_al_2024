List of scripts
* 0_assignment.py (read assignment to open reading frames an RPKM / RPM calculation)
* 1_CI_calc.py (confidence interval calculation)
* 2_onset.py (calculation of chaperone onsets)
* 3_metagene.py (calculation and plotting of metagene profiles of chaperone enrichments)
* 4_cluster.py (DTW based hierarchical substrate clustering)
* 5_domains.py (obtaining CATH annotated domain from the FunFHMMer web server, based on the gene amino acid sequence)
* 6_folding_model.py (alphafold data extraction and calculation of chaperone binding predictions assuming nascent chain folding)
* 7_molten_globule_model.py (alphafold data extraction and calculation of chaperone binding predictions assuming formation of a molten globule)

Reference files
* E. coli MC4100 CDS reference (NZ_HG738867.1.gff3)
* Alphafold files, https://alphafold.ebi.ac.uk/download (UP000000625_83333_ECOLI_v3)
* E. coli MC4100 amino acid sequences
* CATH files are obtained using the FunFHMMer web server
* PAE files are obtained via web scraping with the python “urllib.request” module

1. System Requirements
Hardware Requirements
Normal desktop computer, at least 8 GB of RAM. For faster cluster analysis, a multi-core processor is helpful.

Software Requirements
Python 3.9 with the following additional modules: numpy, statsmodels, scipy, matplotlib, tslearn

2. Installation Guide
Via conda:
Install your conda environment following: https://conda.io/projects/conda/en/latest/user-guide/getting-started.html(https://conda.io/projects/conda/en/latest/user-guide/getting-started.html)
Install python and all modules, for example scipy following: https://anaconda.org/anaconda/scipy(https://anaconda.org/anaconda/scipy)
Installation usually take 2-3 hours. 

3. Demo
Instructions
The demo data set are genome aligned .sam files (the typical output of the Bowtie1 aligner) of an IP (IP_demo.sam) and one total translatome (total_demo.sam). Both files contain 100.000 reads. Scripts are numberd 0-6, corresponding to the respective output demo files. Note that the cluster file output is usually > 100MB and can not be stored in this repository

Expected output
See files with the corresponding number in the demo folder

Expected runtime
All scripts run < 1 minute with the demo dataset except the „similarity_matrix_maker“ function in 4_cluster.py, which usually takes 1 - 2 h with 4000 genes. Run times exponentially increase with genome size, which makes parallel processing favorable with larger genomes. 

4. Instructions for use
How to run the software on your data
Activate the conda environment where python and all modules are installed. To test all pythonscripts, the “path” or “main_path” variable should be a string with the path to the demo folder.

Reproductive instructions
CATH and Alphafold library queries might change over time due to newly deposited data sets. 


