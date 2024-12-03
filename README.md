# Topics Seminar Fall 2024
## https://github.com/EEsdaile/Topics-Seminar-Fall-2024
### This is my final project for the MIP Topics class

The goal is to create a VCF from publicly available BAMs from ancient and modern horses.\
I am attempting to follow the methods from Librado et al. 2024.

Librado, P., Tressières, G., Chauvey, L. et al. Widespread horse-based mobility arose around 2200 bce in Eurasia. Nature 631, 819–825 (2024). https://doi.org/10.1038/s41586-024-07597-5

### <ins>Set up:</ins>
 To set my command line as I like it \
GNU bash, version 4.4.20(1)-release (x86_64-redhat-linux-gnu) \
Copyright (C) 2016 Free Software Foundation, Inc. \
License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>

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

### <ins>Data Download</ins>
Data was downloaded from the European Nucleotide Archive (ENA) under project number PRJEB44430 (https://www.ebi.ac.uk/ena/browser/view/PRJEB44430?show=analyses) using the ENA generated script to download all files. 

This is a long file with a command for each of the 500+ files being downloaded.

GNU Wget 1.19.5 built on linux-gnu.
``` bash
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/{sample_ID}/{file_name}
```
See related files in "download_bams" directory here on GitHub


\
Anaconda was used to create a conda environment with NCBI Datasets to enable downloading of the reference genome (EquCab3). 

anaconda Command line client (version 1.12.1) \
datasets version: 16.33.0 
``` bash    
conda create -n ncbi_datasets  #Create conda environment
conda activate ncbi_datasets   #Activate the new conda environment
conda install -c conda-forge ncbi-datasets-cli  #install NCBI Datasets


#Download the reference genome using NCBI Datasets
datasets download genome accession GCF_002863925.1 --include gff3,rna,cds,protein,genome,seq-report   
```

I indexed the reference genome fasta file using samtools faidx. \
samtools 1.16.1 \
Using htslib 1.16 \
Copyright (C) 2022 Genome Research Ltd. 

``` bash
module load samtools #this enables samtools to be used on the cluster
samtools faidx GCF_002863925.1_EquCab3.0_genomic.fna  #actual command to index the fasta (in this case, an "fna" file)
```

### <ins>Variant Calling</ins> (make_vcf_glf.sh)
ANGSD to was used call SNPs : https://github.com/ANGSD/angsd \
angsd version: 0.941-22-gc877e7f (htslib: 1.21-14-gcf0e7568) build(Oct 30 2024 12:15:14)


ANGSD was used because the quality and coverage of sequencing data from ancient samples, hundreds to thousands of years old, is very poor compared to modern, high quality samples. As such, ANGSD calculates genotypic likelihood and imputation to identify variants across a genome. This does complicate future analysis because it does not produce vcf files in the typical format (0/0, 0/1, 1/1, ./., etc) and instead calculates the likelihood that a sample has each of the possible genotypes (0/0, 0/1, or 1/1). However, there are some custom tools to still conduct common analyses using ancient samples. 


-SNP_pval 0 was used by Librado et al. 2024, but was removed from the code due to it 
causing the script to fail. I used an updated version of ANGSD, which is not the same version 
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


### <ins>Output</ins>
The output file was unzipped to view the contents \
gunzip (gzip) 1.9 \
Copyright (C) 2007, 2011-2017 Free Software Foundation, Inc.

head (GNU coreutils) 8.30 \
Copyright (C) 2018 Free Software Foundation, Inc.
``` bash
gunzip -k output2.bcf.beagle.gz #unzip file and keep original file (-k)
head -n 20 output2.bcf.beagle # view first 20 lines of the unzipped file.
```
This is the output of the head of the file. 
![image](https://github.com/user-attachments/assets/1543fabe-f80e-4c5c-aa42-9e94c65a80c2)

---------------------------------------------------------------------------------

### <ins> What next? </ins>
The beagle output file can be fed into PCAngsd to generate a PCA to evaluate how ancient and modern horse and other equid samples cluster together. For example, I would expect the modern european horses to cluster with the ancient european horses, and the ancient eastern asian horses to cluster with modern eastern asian horses. Additional analyses could be done to evaluate gene flow. Alternatively, modern horses may have drifted away from ancient horses. Regardless, it is likely that the clustering would show and interesting pattern due to human selection over centuries. 
https://www.popgen.dk/software/index.php/PCAngsd

This would be a starting place to run the PCA: \
PCAngsd **Version 1.35**
``` bash
#!/bin/bash

#SBATCH --partition=amilan
#SBATCH --nodes=2
#SBATCH --output=pcangsd
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=e.s.esdaile@colostate.edu


pcangsd -beagle output2.bcf.beagle.gz -e 20 -t 64 -o /scratch/alpine/c829993361@colostate.edu/ancient_horses/PCA_output_20241125
```

### <ins>Reproducibility</ins>
Key steps to include reproducibility include:
1. Persistent data sharing - Publically available and accessible data. Ideally with properly logged metadata so that the data is findable without having to search through every publication. The data I used was already publicly available in a repository commonly used by the equine genetics community.
2. Code version control and persistent sharing -  Including the version of all software used and making the code publicly available in a repository. I believe that I have done so by uploading my code to GitHub.
3. Literate programming - describe what each step does and don't use mysterious variables. I have tried to do this because I know if I revisit this in the future, I won't remember what any of the steps are, so how can I expect someone else to follow them for the first time.
4. Documentation - This README :)
5. Compute environment and control - I didn't do as good of a job with this, but do have one step in a conda environment. Ideally I would load it all into a conda environment, but some programs were downloaded from github and I need to learn how to integrate those into a conda environment.

### <ins> Peer Feedback </ins>
I didn't receive much feedback, and the little I received was related to version control and literate programming. Usually my code is a bit of a mess as I write it and work through the problems and errors, and once I have it sorted out, I go back through and clean everything up, ensuring that it still works along the way. I have since done that and am happy with my documentation and code, and have included the versions of all the softwares used. 

