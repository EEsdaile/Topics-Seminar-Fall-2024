#!/bin/bash

#SBATCH --partition=amilan
#SBATCH --nodes=6
#SBATCH --output=download-PRJEB44430
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=e.s.esdaile@colostate.edu

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3038933/Gral10_AMIS-1-01745_Fra_m18924.Horse_mtDNA_30firstbp_copiedattheend.realigned.q25.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3034591/RusxNx113_AMIS-1-00912_Rus_m13978.Bowtie2.EquCab3.pmds1.r.t.m.s.q25.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3034540/NB44_CGG-1-017036_Rus_m1856.Bowtie2.EquCab3.pmds1.r.t.m.s.q25.bam.bai
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3039045/RusxNx129_AMIS-1-00928_Rus_m36541.Horse_mtDNA_30firstbp_copiedattheend.realigned.q25.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3038952/Hasanlu1140_CGG-1-019998_Ira_m663.Horse_mtDNA_30firstbp_copiedattheend.realigned.q25.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3034579/Rus2017x5_AMIS-1-00776_Rus_m51700i.Bowtie2.EquCab3.pmds1.r.t.m.s.q25.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3038940/Gral5_AMIS-1-01740_Fra_m18924.Horse_mtDNA_30firstbp_copiedattheend.realigned.q25.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3034564/PLMoz1x1_AMIS-1-01703_Pol_m3801.Bowtie2.EquCab3.pmds1.r.t.m.s.q25.bam.bai
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3034429/EREsto3_CGG-1-021522_Est_m659.Bowtie2.EquCab3.pmds1.r.t.m.s.q25.bam.bai
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3039033/RusxNx02_AMIS-1-00564_Rus_m1825.Horse_mtDNA_30firstbp_copiedattheend.realigned.q25.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3034550/Novoil2_Ages-1-7476_Kaz_m1832.Bowtie2.EquCab3.pmds1.r.t.m.s.q25.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3034617/Tarpan_AMIS_1_02430_Ukr_1868.Bowtie2.EquCab3.pmds1.r.t.m.s.q25.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3038938/Gral3_AMIS-1-01738_Fra_m18924.Horse_mtDNA_30firstbp_copiedattheend.realigned.q25.bam.bai
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3039068/Ural2017x29_AMIS-1-00418_Rus_m1853.Horse_mtDNA_30firstbp_copiedattheend.realigned.q25.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3034601/RusxNx53_AMIS-1-00616_Rus_m2785.Bowtie2.EquCab3.pmds1.r.t.m.s.q25.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3038987/LORUS2018x43_AMIS-1-01311_Rus_m1798.Horse_mtDNA_30firstbp_copiedattheend.realigned.q25.bam.bai
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3034603/RusxNx96_AMIS-1-00895_Rus_m4586.Bowtie2.EquCab3.pmds1.r.t.m.s.q25.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3034617/Tarpan_AMIS_1_02430_Ukr_1868.Bowtie2.EquCab3.pmds1.r.t.m.s.q25.bam.bai
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3034457/GVA9035_CGG-1-021614_Mon_m1024.Bowtie2.EquCab3.pmds1.r.t.m.s.q25.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3038900/Botai2018x26_AMIS-1-01642_Kaz_m3333.Horse_mtDNA_30firstbp_copiedattheend.realigned.q25.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3034633/Ural2017x29_AMIS-1-00418_Rus_m1853.Bowtie2.EquCab3.pmds1.r.t.m.s.q25.bam.bai
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3039041/RusxNx114_AMIS-1-00913_Rus_m15249.Horse_mtDNA_30firstbp_copiedattheend.realigned.q25.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3038953/Hasanlu2327_CGG-1-019995_Ira_m768.Horse_mtDNA_30firstbp_copiedattheend.realigned.q25.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3034582/RusxNx02_AMIS-1-00564_Rus_m1825.Bowtie2.EquCab3.pmds1.r.t.m.s.q25.bam.bai
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3038991/Mon2017x126x1_CGG-1-021881_Mon_m942.Horse_mtDNA_30firstbp_copiedattheend.realigned.q25.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3034536/Monx9_AMIS-1-01188_Mon_m775.Bowtie2.EquCab3.pmds1.r.t.m.s.q25.bam.bai
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3034568/PRAGUE32_AMIS-1-00842_Cze_m787.Bowtie2.EquCab3.pmds1.r.t.m.s.q25.bam.bai
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3034475/Hohler1x2_AMIS_1_02432_Ger_CWC.Bowtie2.EquCab3.pmds1.r.t.m.s.q25.bam.bai
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3038923/DJM613x1_AMIS-1-02451_Dan_m2994.Horse_mtDNA_30firstbp_copiedattheend.realigned.q25.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3038981/LORUS2018x09_AMIS-1-01277_Rus_m3131.Horse_mtDNA_30firstbp_copiedattheend.realigned.q25.bam.bai
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3038925/EREsto2_CGG-1-021521_Est_m661.Horse_mtDNA_30firstbp_copiedattheend.realigned.q25.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3034446/Gral1_AMIS-1-01736_Fra_m16101.Bowtie2.EquCab3.pmds1.r.t.m.s.q25.bam.bai
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3034499/KZ2019x1a_AMIS-1-02014_Kaz_m293.Bowtie2.EquCab3.pmds1.r.t.m.s.q25.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3034657/Var3_AMIS-1-01065_Rus_m5549.Bowtie2.EquCab3.pmds1.r.t.m.s.q25.bam.bai
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3038913/BPTDG1_CGG-1-020427_Fra_UPal.Horse_mtDNA_30firstbp_copiedattheend.realigned.q25.bam.bai
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3038815/Ural2017x20_AMIS-1-00409_Rus_m1356.Horse_mtDNA_30firstbp_copiedattheend.realigned.q25.bam.bai
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3034546/Nov6x6_AMIS-1-02629_Rus_m3224.Bowtie2.EquCab3.pmds1.r.t.m.s.q25.bam.bai
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3034495/KSH5_CGG-1-017099_Kaz_m1895.Bowtie2.EquCab3.pmds1.r.t.m.s.q25.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3034378/Bi5x1_AMIS-1-02594_Rus_m625.Bowtie2.EquCab3.pmds1.r.t.m.s.q25.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3034639/Ural2017x36_AMIS-1-00425_Rus_m1892.Bowtie2.EquCab3.pmds1.r.t.m.s.q25.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3039067/Ural2017x26_AMIS-1-00415_Rus_m2751.Horse_mtDNA_30firstbp_copiedattheend.realigned.q25.bam.bai
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3038969/KB221_AMIS-1-01629_Rus_PAL.Horse_mtDNA_30firstbp_copiedattheend.realigned.q25.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3038950/Halvai2_AMIS-1-00058_Kaz_m1856.Horse_mtDNA_30firstbp_copiedattheend.realigned.q25.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3034418/Cho1x3_AMIS-1-02597_Rus_m625.Bowtie2.EquCab3.pmds1.r.t.m.s.q25.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3038784/LORUS2018x19_AMIS-1-01287_Rus_m5391.Horse_mtDNA_30firstbp_copiedattheend.realigned.q25.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3034382/Borly8_AMIS-1-00040_Kaz_m4290.Bowtie2.EquCab3.pmds1.r.t.m.s.q25.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3034498/KY001_CGG-1-022896_Kaz_m3244.Bowtie2.EquCab3.pmds1.r.t.m.s.q25.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3034388/Botai2018x15_AMIS-1-01653_Kaz_m3333.Bowtie2.EquCab3.pmds1.r.t.m.s.q25.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3038761/AC7970_AMIS-1-00131_Tur_m290.Horse_mtDNA_30firstbp_copiedattheend.realigned.q25.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3038911/BotaixI_CGG-1-020189_Kaz_m3333.Horse_mtDNA_30firstbp_copiedattheend.realigned.q25.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3034480/Hohler3x3_AMIS_1_02436_Ger_CWC.Bowtie2.EquCab3.pmds1.r.t.m.s.q25.bam.bai
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3038903/Botai2018x31_AMIS-1-01659_Kaz_m3333.Horse_mtDNA_30firstbp_copiedattheend.realigned.q25.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3034417/Chin1x1_AMIS-1-01232_Rus_m810.Bowtie2.EquCab3.pmds1.r.t.m.s.q25.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3038766/Bi5x1_AMIS-1-02594_Rus_m625.Horse_mtDNA_30firstbp_copiedattheend.realigned.q25.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3038971/Kent6_AMIS-1-0149_Kaz_m1463.Horse_mtDNA_30firstbp_copiedattheend.realigned.q25.bam.bai
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3034440/Gar3_CGG-1-018389_Rom_m1539.Bowtie2.EquCab3.pmds1.r.t.m.s.q25.bam.bai
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3034496/KSK16232_AMIS-1-00027_Tur_m6319.Bowtie2.EquCab3.pmds1.r.t.m.s.q25.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3038974/KSH5_CGG-1-017099_Kaz_m1895.Horse_mtDNA_30firstbp_copiedattheend.realigned.q25.bam.bai
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3039052/ROxCASx10_CGG-1-022929_Rom_m4309.Horse_mtDNA_30firstbp_copiedattheend.realigned.q25.bam.bai
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3038812/Tarpan_AMIS_1_02430_Ukr_1868.Horse_mtDNA_30firstbp_copiedattheend.realigned.q25.bam.bai
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3034649/Ural2017x69_AMIS-1-00458_Rus_m2774.Bowtie2.EquCab3.pmds1.r.t.m.s.q25.bam.bai
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3038794/Monx7_AMIS-1-01186_Mon_m1030.Horse_mtDNA_30firstbp_copiedattheend.realigned.q25.bam.bai
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3034446/Gral1_AMIS-1-01736_Fra_m16101.Bowtie2.EquCab3.pmds1.r.t.m.s.q25.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3034570/PRAGUE40_AMIS-1-00850_Cze_m2037.Bowtie2.EquCab3.pmds1.r.t.m.s.q25.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3034524/Mon2017x176_CGG-1-0219435_Mon_m939.Bowtie2.EquCab3.pmds1.r.t.m.s.q25.bam.bai
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3034524/Mon2017x176_CGG-1-0219435_Mon_m939.Bowtie2.EquCab3.pmds1.r.t.m.s.q25.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3034578/Rus2017x3_CGG-NO-NUMBER_Rus_m51700i.Bowtie2.EquCab3.pmds1.r.t.m.s.q25.bam.bai
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3034535/Monx7_AMIS-1-01186_Mon_m1030.Bowtie2.EquCab3.pmds1.r.t.m.s.q25.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3038801/PLKaz4_AMIS-1-01680_Pol_m1550.Horse_mtDNA_30firstbp_copiedattheend.realigned.q25.bam.bai
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3034449/Gral4_AMIS-1-01739_Fra_m18924.Bowtie2.EquCab3.pmds1.r.t.m.s.q25.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3039029/Rus2017x3_CGG-NO-NUMBER_Rus_m51700i.Horse_mtDNA_30firstbp_copiedattheend.realigned.q25.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3034608/ERRus30x31_CGG-1-021550x1_Rus_m3435.Bowtie2.EquCab3.pmds1.r.t.m.s.q25.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3038973/Kozhai_AMIS-1-00155x051_Kaz_m3342.Horse_mtDNA_30firstbp_copiedattheend.realigned.q25.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3039040/RusxNx113_AMIS-1-00912_Rus_m13978.Horse_mtDNA_30firstbp_copiedattheend.realigned.q25.bam.bai
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3038887/Besta7_AMIS-1-00063_Kaz_m3132.Horse_mtDNA_30firstbp_copiedattheend.realigned.q25.bam.bai
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3034397/Botai2018x3_AMIS-1-01643_Kaz_m3333.Bowtie2.EquCab3.pmds1.r.t.m.s.q25.bam.bai
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3034591/RusxNx113_AMIS-1-00912_Rus_m13978.Bowtie2.EquCab3.pmds1.r.t.m.s.q25.bam.bai
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3034603/RusxNx96_AMIS-1-00895_Rus_m4586.Bowtie2.EquCab3.pmds1.r.t.m.s.q25.bam.bai
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3038919/Chalk3_AMIS-1-01390_UK_m474.Horse_mtDNA_30firstbp_copiedattheend.realigned.q25.bam.bai
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3038782/Kent2_AMIS-1-00145_Kaz_m1474.Horse_mtDNA_30firstbp_copiedattheend.realigned.q25.bam.bai
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3038766/Bi5x1_AMIS-1-02594_Rus_m625.Horse_mtDNA_30firstbp_copiedattheend.realigned.q25.bam.bai
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3034509/LORUS2018x23_AMIS-1-01291_Rus_m1796.Bowtie2.EquCab3.pmds1.r.t.m.s.q25.bam.bai
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3039038/RusxNx111_AMIS-1-00910_Rus_m22013.Horse_mtDNA_30firstbp_copiedattheend.realigned.q25.bam.bai
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3034372/Arz1x12_AMIS-1-01222_Rus_m854.Bowtie2.EquCab3.pmds1.r.t.m.s.q25.bam.bai
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3038756/PAVH8_CGG-1-018165_Kaz_m2961.Bowtie2.EquCab3.pmds1.r.t.m.s.q25.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3034628/Ural2017x20_AMIS-1-00409_Rus_m1356.Bowtie2.EquCab3.pmds1.r.t.m.s.q25.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3034659/Zam9_CGG-1-018369_Por_m2559.Bowtie2.EquCab3.pmds1.r.t.m.s.q25.bam.bai
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3034473/HasanluV31E_CGG-1-021461_Ira_m768.Bowtie2.EquCab3.pmds1.r.t.m.s.q25.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3038764/Arz1x13_AMIS-1-01223_Rus_m581.Horse_mtDNA_30firstbp_copiedattheend.realigned.q25.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3038776/Ganx10_AMIS-1-01210_Mon_m775.Horse_mtDNA_30firstbp_copiedattheend.realigned.q25.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3034399/Botai2018x6_AMIS-1-01646_Kaz_m3333.Bowtie2.EquCab3.pmds1.r.t.m.s.q25.bam.bai
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3038809/Sho1x2_AMIS-1-02612_Rus_m7500.Horse_mtDNA_30firstbp_copiedattheend.realigned.q25.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3038910/BotaixG_CGG-1-020187_Kaz_m3193.Horse_mtDNA_30firstbp_copiedattheend.realigned.q25.bam.bai
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3034574/PRAGUE77_AMIS-1-01716_Cze_m624.Bowtie2.EquCab3.pmds1.r.t.m.s.q25.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3038754/SAGxS27_CGG-1-019559_Ira_m1102.Bowtie2.EquCab3.pmds1.r.t.m.s.q25.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3038805/PRAGUE40_AMIS-1-00850_Cze_m2037.Horse_mtDNA_30firstbp_copiedattheend.realigned.q25.bam.bai
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3034598/RusxNx129_AMIS-1-00928_Rus_m36541.Bowtie2.EquCab3.pmds1.r.t.m.s.q25.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3034418/Cho1x3_AMIS-1-02597_Rus_m625.Bowtie2.EquCab3.pmds1.r.t.m.s.q25.bam.bai
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3034532/MOLDA1_CGG-1-022639_Mol_m2063.Bowtie2.EquCab3.pmds1.r.t.m.s.q25.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3039065/Ural2017x1_AMIS-1-00391_Rus_m1868.Horse_mtDNA_30firstbp_copiedattheend.realigned.q25.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3034636/Ural2017x31_CGG-NO-NUMBER_Rus_m1845.Bowtie2.EquCab3.pmds1.r.t.m.s.q25.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3039034/RusxNx03_AMIS-1-00565_Rus_m1851.Horse_mtDNA_30firstbp_copiedattheend.realigned.q25.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3034403/BotaixK_CGG-1-020191_Kaz_m3501.Bowtie2.EquCab3.pmds1.r.t.m.s.q25.bam.bai
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3034484/KB218_AMIS-1-01627_Rus_m42475.Bowtie2.EquCab3.pmds1.r.t.m.s.q25.bam.bai
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3034366/18ELTu18_AMIS-1-01102_Spa_m588.Bowtie2.EquCab3.pmds1.r.t.m.s.q25.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3038974/KSH5_CGG-1-017099_Kaz_m1895.Horse_mtDNA_30firstbp_copiedattheend.realigned.q25.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3038982/LORUS2018x12_AMIS-1-01280_Rus_m2995.Horse_mtDNA_30firstbp_copiedattheend.realigned.q25.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3039083/Ural2017x80_AMIS-1-01407_Kaz_m3235.Horse_mtDNA_30firstbp_copiedattheend.realigned.q25.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3039069/Ural2017x30_AMIS-1-00419_Rus_m1863.Horse_mtDNA_30firstbp_copiedattheend.realigned.q25.bam.bai
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3034503/LORUS2018x09_AMIS-1-01277_Rus_m3131.Bowtie2.EquCab3.pmds1.r.t.m.s.q25.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3038813/Tarquinia3298_AMIS-1-01381_Ita_m657.Horse_mtDNA_30firstbp_copiedattheend.realigned.q25.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ303/ERZ3034637/Ural2017x33_AMIS-1-00422_Rus_m1854.Bowtie2.EquCab3.pmds1.r.t.m.s.q25.bam.bai
