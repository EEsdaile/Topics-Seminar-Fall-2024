# Topics Seminar Fall 2024
## This is my final project for the MIP Topics class

The goal is to create a VCF from publicly available BAMs from ancient and modern horses. 
From this VCF I will create a structure plot and a PCA by breed to see if the ancient horses cluster as expected. 
I am attempting to follow the methods from the ancient horse papers. 

 To set my command line as I like it
 ` bash
 PS1='\u:\w\$ '
 `
 
 To set sq to print my squeue 
 ` bash
 function sq() {
    local user=${1:-$(whoami)}
    cmd="squeue -u ${user}"
    echo "# ${cmd}"
    $cmd
}
`

Paper: https://www.nature.com/articles/s41586-024-07597-5#MOESM1
Data was downloaded from https://www.ebi.ac.uk/ena/browser/view/PRJEB44430?show=analyses using the providing script to download all files. 

EquCab3 downloaded with NCBI Datasets
	In conda environment
` bash
conda create -n ncbi_datasets
conda activate ncbi_datasets
conda install -c conda-forge ncbi-datasets-cli

datasets download genome accession GCF_002863925.1 --include gff3,rna,cds,protein,genome,seq-report
`

ANGSD to call SNPs : https://github.com/ANGSD/angsd
-SNP_pval 0 was used per the original publication


