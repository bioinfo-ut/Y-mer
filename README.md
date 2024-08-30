# Y-mer
Predicting human Y chromosome haplogroups from ultra-low depth sequencing data 


#COUNTING K-MER FREQUENCIES
#GenomeTester4 package contains programs gmer_counter, glistmaker, glistquery https://github.com/bioinfo-ut/Genometester4


With fastq or fasta files:

gmer_counter -dbb k-mers.dbb sample.fastq |cut -f 3 |tail -n +3 > sample_G.counts


With bam or cram files:

samtools fasta sample.bam|gmer_counter -dbb k-mers.dbb - |cut -f 3 |tail -n +3 > sample_G.counts


With having GenomeTester4 based list file, mandatory if using multiple models:

glistquery sample_25.list -f k-mers.txt |cut -f 2 > sample_G.counts



#CALLING HG-s
 
Calling commandline order is script, model file, sample counts file and R formated output file name

Rscript PREDICTER.R model.Rdata sample_G.counts sample_result.RData > sample_result.txt


#SAMPLES

aDNA sample DA189 fastq reads ERR2505887.fastq.gz
aDNA sample DA189 mapped bam file DA189.sort.rmdup.realign.md.bam

assembled chrY NA20509.HIFIRW.ONTUL.na.chrY.fasta

#MODELS

M213E
https://bioinfo.ut.ee/randomtandem/mudelid/

modelfile                                M213E/M213E_50k.Rdata
k-mer dictionary for glistquery          M213E/M213E_50k.txt
k-mer binary dictionary for gmer_counter M213E/M213E_50k.dbb

#WEB tool

M21E
https://bioinfo.ut.ee/randomtandem/magic/
