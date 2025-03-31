**# Y-mer**

**Y-mer is a workflow for determining the Y chromosome haplogroup from very low-coverage sequencing data. 
For this purpose, we have developed a workflow that allows users to create their own model from 
high-coverage sequencing data or use models provided by us. 
Our models have a sensitivity on 0.005x sequencing coverage in the Y chromosome and even lower, depending on the haplogroup**

Please cite: Puurand, T. et al. (2024) ‘Y-mer: A k-mer based method for determining human Y chromosome haplogroups from 
ultra-low sequencing depth data’. Available at: https://doi.org/10.21203/rs.3.rs-5042960/v1

Required programs:\
R\
perl\
GenomeTester4 https://github.com/bioinfo-ut/Genometester4 

Please download files from directories 'model_training' or 'haplogroup_prediction' and https://doi.org/10.5281/zenodo.15089783
according your interest.


**MODEL TRAINING**

Insert samples data in files men.txt and women.txt, create directories 'data' and 'lists', look over 'bam' files location ans run
'**perl Y-mer.pl**'. SDD requirements is 30 GB per sample and RAM 80 GB during processing.


**HAPLOGROUP PREDICTION**

A. COUNTING K-MER FREQUENCIES

With fastq or fasta files:

gmer_counter -dbb k-mers.dbb sample.fastq |cut -f 3 |tail -n +3 > sample_G.counts


With bam or cram files:

samtools fasta sample.bam|gmer_counter -dbb k-mers.dbb - |cut -f 3 |tail -n +3 > sample_G.counts


With having GenomeTester4 based list file, mandatory if using multiple models:

**glistquery sample_25.list -f k-mers.txt |cut -f 2 > sample_G.counts**



B. CALLING HG-s
 
Calling commandline order is script model file sample counts file and R formated output file name

**Rscript PREDICTER.R model.Rdata sample_G.counts sample_result.RData > sample_result.txt**


#WEB tool
https://bioinfo.ut.ee/randomtandem/Y-mer/


#SAMPLES

aDNA sample DA189 fastq reads https://bioinfo.ut.ee/randomtandem/mudelid/ERR2505887.fastq.gz
aDNA sample DA189 mapped bam file https://bioinfo.ut.ee/randomtandem/mudelid/DA189.sort.rmdup.realign.md.bam

assembled chrY https://bioinfo.ut.ee/randomtandem/mudelid/NA20509.HIFIRW.ONTUL.na.chrY.fasta



#MODELS

https://bioinfo.ut.ee/randomtandem/mudelid/ and https://doi.org/10.5281/zenodo.15089783

modelfile                                M213E/M213E_50k.Rdata
k-mer dictionary for glistquery          M213E/M213E_50k.txt
k-mer binary dictionary for gmer_counter M213E/M213E_50k.dbb
