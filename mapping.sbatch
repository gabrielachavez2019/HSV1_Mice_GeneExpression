#!/bin/bash
#SBATCH -p batch
#SBATCH -t 5-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mail-user=gabriela.toomer@okstate.edu
#SBATCH --mail-type=end
module load bowtie2/2.3.4.1 fastqc/0.11.7 samtools/1.10 trimmomatic/0.38 cufflinks/2.2.1 tophat/2.1.1 R/4.0.3

##Concatenate files Mocks or uninfected:
#  echo Begining concatenating files...
#  
#cat /scratch/gatoo/SecNGS/LIB109080_S*_R1_001.fastq.gz >> /scratch/gatoo/HSV1-mouse/Fem-dLAT_F.fastq.gz
#cat /scratch/gatoo/SecNGS/LIB109080_S*_R2_001.fastq.gz >> /scratch/gatoo/HSV1-mouse/Fem-dLAT_R.fastq.gz
#
#cat /scratch/gatoo/SecNGS/LIB109081_S*_R1_001.fastq.gz >> /scratch/gatoo/HSV1-mouse/Fem-WT_F.fastq.gz
#cat /scratch/gatoo/SecNGS/LIB109081_S*_R2_001.fastq.gz >> /scratch/gatoo/HSV1-mouse/Fem-WT_R.fastq.gz
#
#cat /scratch/gatoo/SecNGS/LIB109082_S*_R1_001.fastq.gz >> /scratch/gatoo/HSV1-mouse/Fem-U1_F.fastq.gz
#cat /scratch/gatoo/SecNGS/LIB109082_S*_R2_001.fastq.gz >> /scratch/gatoo/HSV1-mouse/Fem-U1_R.fastq.gz
#
#cat /scratch/gatoo/SecNGS/LIB109083_S*_R1_001.fastq.gz >> /scratch/gatoo/HSV1-mouse/Male-U1_F.fastq.gz
#cat /scratch/gatoo/SecNGS/LIB109083_S*_R2_001.fastq.gz >> /scratch/gatoo/HSV1-mouse/Male-U1_R.fastq.gz
#
#cat /scratch/gatoo/SecNGS/LIB109084_S*_R1_001.fastq.gz >> /scratch/gatoo/HSV1-mouse/Male-dLAT_F.fastq.gz
#cat /scratch/gatoo/SecNGS/LIB109084_S*_R2_001.fastq.gz >> /scratch/gatoo/HSV1-mouse/Male-dLAT_R.fastq.gz
#
#cat /scratch/gatoo/SecNGS/LIB109085_S*_R1_001.fastq.gz >> /scratch/gatoo/HSV1-mouse/Male-WT_F.fastq.gz
#cat /scratch/gatoo/SecNGS/LIB109085_S*_R2_001.fastq.gz >> /scratch/gatoo/HSV1-mouse/Male-WT_R.fastq.gz

#BoHSV1-cow
#For this part you need to be in /home/gatoo/HSV1-mouse
#echo Finishing full size library files...
##echo Bulding HSV1 genome ...
#bowtie2-build --num-threads 32 /scratch/gatoo/HSV1-mouse/HSV1.fasta  /scratch/gatoo/HSV1-mouse/HSV1

names='Fem-dLAT Fem-WT Fem-U1
       Male-dLAT Male-WT Male-U1'

for name in $names
 do
  #echo Starting mapping of "$name" against the virus
  echo Starting mapping of "$name" against the mouse genome
#tophat2 --num-threads 32 -o /scratch/gatoo/HSV1-mouse/output_"$name"  -G /scratch/gatoo/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf --transcriptome-index /scratch/gatoo/Mus_musculus/UCSC/mm10/Annotation/Genes/ --read-edit-dist 2 --read-gap-length 2 --read-mismatches 2 --read-realign-edit-dist 1 --library-type fr-unstranded --min-intron-length 50 --max-intron-length 300000 /scratch/gatoo/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome /scratch/gatoo/HSV1-mouse/"$name"_F.fastq.gz /scratch/gatoo/HSV1-mouse/"$name"_R.fastq.gz

tophat2 --num-threads 32 -o /scratch/gatoo/HSV1-mouse/output_HSV1_"$name" /scratch/gatoo/HSV1-mouse/HSV1 /scratch/gatoo/HSV1-mouse/"$name"_F.fastq.gz /scratch/gatoo/HSV1-mouse/"$name"_R.fastq.gz


cufflinks --num-threads 32 -o /scratch/gatoo/HSV1-mouse/output_cuff_HSV1_"$name"  --frag-bias-correct /scratch/gatoo/HSV1-mouse/HSV1.fasta --multi-read-correct -G /scratch/gatoo/HSV1-mouse/HSV1.gff3 /scratch/gatoo/HSV1-mouse/output_HSV1_"$name"/accepted_hits.bam

cp  /scratch/gatoo/HSV1-mouse/output_cuff_HSV1_"$name"/genes.fpkm_tracking HSV1_"$name".txt

   echo $name counting: done!
  done
echo All done
