#!/bin/bash

#SBATCH --partition=amilan
#SBATCH --nodes=1
#SBATCH --output=angsd_make_vcf_glf
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=e.s.esdaile@colostate.edu


#call variants with angsd
#-snp_pval 0 per original publicaion

#/projects/c829993361@colostate.edu/software/angsd/angsd -b bams_full_path_test10.txt\
# -dovcf 1 -gl 1 -dopost 1 -domajorminor 1 -domaf 1\
# -minQ 30 -minMapQ 25 -uniqueOnly 1 -remove_bads 1 -C 50 -snp_pval 0\
# -ref /projects/c829993361@colostate.edu/EquCab3/ncbi_dataset/data/GCF_002863925.1/GCF_002863925.1_EquCab3.0_genomic.fna\
# -fai /projects/c829993361@colostate.edu/EquCab3/ncbi_dataset/data/GCF_002863925.1/GCF_002863925.1_EquCab3.0_genomic.fna.fai

/projects/c829993361@colostate.edu/software/angsd/angsd -b bams_full_path_test10.txt\
 -dobcf 1 -gl 1 -doGlf 2 -dopost 1 -domajorminor 1 -domaf 1\
 -minQ 30 -minmapQ 25 -uniqueonly 1 -remove_bads 1 -C 50\
 -out output2.bcf \
 -ref /projects/c829993361@colostate.edu/EquCab3/ncbi_dataset/data/GCF_002863925.1/GCF_002863925.1_EquCab3.0_genomic_Chr_names.fna\
 -fai /projects/c829993361@colostate.edu/EquCab3/ncbi_dataset/data/GCF_002863925.1/GCF_002863925.1_EquCab3.0_genomic_Chr_names.fna.fai

