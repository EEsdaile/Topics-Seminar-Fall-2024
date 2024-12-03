# Topics Seminar Fall 2024
## This is my final project for the MIP Topics class

The goal is to create a VCF from publicly available BAMs from ancient and modern horses.\
I am attempting to follow the methods from Librado et al. 2024.

Librado, P., Tressières, G., Chauvey, L. et al. Widespread horse-based mobility arose around 2200 bce in Eurasia. Nature 631, 819–825 (2024). https://doi.org/10.1038/s41586-024-07597-5

 To set my command line as I like it 
 ``` bash  
 PS1='\u:\w\$ '  
 ```
 
 To set sq as a shortcut to print my squeue
 ``` bash   
 function sq() {   
    local user=${1:-$(whoami)}   
    cmd="squeue -u ${user}"   
    echo "# ${cmd}"   
    $cmd   
}   
```

Data was downloaded from the European Nucleotide Archive (ENA) under project number PRJEB44430 (https://www.ebi.ac.uk/ena/browser/view/PRJEB44430?show=analyses) using the ENA generated script to download all files.

``` bash
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/{sample_ID}/{file_name}
```

The reference genome, EquCab3, downloaded with NCBI Datasets \
using a conda environment 
``` bash    
conda create -n ncbi_datasets  
conda activate ncbi_datasets  
conda install -c conda-forge ncbi-datasets-cli  

datasets download genome accession GCF_002863925.1 --include gff3,rna,cds,protein,genome,seq-report   
```

I indexed the reference genome fasta file using samtools faidx
``` bash
module load samtools
samtools faidx GCF_002863925.1_EquCab3.0_genomic.fna
```

ANGSD to was used call SNPs : https://github.com/ANGSD/angsd   \
-SNP_pval 0 was used per the original publication but was removed from the code due to it \
causing the script to fail. I used an updated version of ANGSD, which is not the same version \
used by Librado et al. 2024. \
The script was run using slurm
``` bash
#!/bin/bash

#SBATCH --partition=amilan
#SBATCH --nodes=1
#SBATCH --output=angsd_make_vcf_glf
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=e.s.esdaile@colostate.edu


#call variants with angsd
#-snp_pval 0 per original publicaion

/projects/c829993361@colostate.edu/software/angsd/angsd -b bams_full_path_test10.txt\
 -dobcf 1 -gl 1 -doGlf 2 -dopost 1 -domajorminor 1 -domaf 1\
 -minQ 30 -minmapQ 25 -uniqueonly 1 -remove_bads 1 -C 50\
 -out output2.bcf \
 -ref /projects/c829993361@colostate.edu/EquCab3/ncbi_dataset/data/GCF_002863925.1/GCF_002863925.1_EquCab3.0_genomic_Chr_names.fna\
 -fai /projects/c829993361@colostate.edu/EquCab3/ncbi_dataset/data/GCF_002863925.1/GCF_002863925.1_EquCab3.0_genomic_Chr_names.fna.fai
```

-b is my file of bams to use (4 horses) \
-dobcf created a bcf file (binary vcf file) \
-gl 1 called genotype likelihoods with model #1 \
-doGlf 2 creats a Beagle output file (calls haplotypes) \
-dopost 1 calculates a posterior probability using the frequency as prior \
-domajorminor 1 infer the major and minor alleles using maximum likelihood. \
-domaf 1 calculates fixed major and minor allele frequencies an outputs a *.maf.gz file \
-minQ 30 uses a minimum base quality (phred score) of 30 \
-minmapQ 25 uses a minimum mapping quality of 25 \
-uniqueonly 1 remove reads with multiple best hits \
-remove_bads 1 removed reads with a flag above 255 \
-C 50 adjusts the mapQ for excessive mismatches \
-ref and -fai are the reference genome and index file for the reference genome \
-out is the prefix for the generated output files. \


