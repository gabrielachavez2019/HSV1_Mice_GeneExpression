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
##
#cat /scratch/gatoo/SecNGS/LIB109081_S*_R1_001.fastq.gz >> /scratch/gatoo/HSV1-mouse/Fem-WT_F.fastq.gz
#cat /scratch/gatoo/SecNGS/LIB109081_S*_R2_001.fastq.gz >> /scratch/gatoo/HSV1-mouse/Fem-WT_R.fastq.gz
##
#cat /scratch/gatoo/SecNGS/LIB109082_S*_R1_001.fastq.gz >> /scratch/gatoo/HSV1-mouse/Fem-U1_F.fastq.gz
#cat /scratch/gatoo/SecNGS/LIB109082_S*_R2_001.fastq.gz >> /scratch/gatoo/HSV1-mouse/Fem-U1_R.fastq.gz
##
#cat /scratch/gatoo/SecNGS/LIB109083_S*_R1_001.fastq.gz >> /scratch/gatoo/HSV1-mouse/Male-U1_F.fastq.gz
#cat /scratch/gatoo/SecNGS/LIB109083_S*_R2_001.fastq.gz >> /scratch/gatoo/HSV1-mouse/Male-U1_R.fastq.gz
##
#cat /scratch/gatoo/SecNGS/LIB109084_S*_R1_001.fastq.gz >> /scratch/gatoo/HSV1-mouse/Male-dLAT_F.fastq.gz
#cat /scratch/gatoo/SecNGS/LIB109084_S*_R2_001.fastq.gz >> /scratch/gatoo/HSV1-mouse/Male-dLAT_R.fastq.gz
##
#cat /scratch/gatoo/SecNGS/LIB109085_S*_R1_001.fastq.gz >> /scratch/gatoo/HSV1-mouse/Male-WT_F.fastq.gz
#cat /scratch/gatoo/SecNGS/LIB109085_S*_R2_001.fastq.gz >> /scratch/gatoo/HSV1-mouse/Male-WT_R.fastq.gz

#echo Finishing building full size library files...

##For this part you need to be in /home/gatoo/HSV1-mouse

#echo Bulding HSV1 genome ...

##Must run from  /scratch/gatoo/HSV1-mouse/ 

#bowtie2-build /scratch/gatoo/HSV1-mouse/HSV1.fasta HSV1

names='Fem-dLAT Fem-WT Fem-U1
       Male-dLAT Male-WT Male-U1'

for name in $names
 do
   echo Read Quality Analysis "$name"
   fastqc "$name"_F.fastq.gz "$name"_R.fastq.gz
   cutadapt "$name"_F.fastq.gz "$name"_R.fastq.gz 
   
   echo Starting mapping of "$name"  against the mouse genome
  
      tophat2 --num-threads 32 -o /scratch/gatoo/HSV1-mouse/output_"$name"  -G /scratch/gatoo/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf --transcriptome-index /scratch/gatoo/Mus_musculus/UCSC/mm10/Annotation/Genes/ --read-edit-dist 2 --read-gap-length 2 --read-mismatches 2 --read-realign-edit-dist 1 --library-type fr-unstranded --min-intron-length 50 --max-intron-length 300000 /scratch/gatoo/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome /scratch/gatoo/HSV1-mouse/"$name"_F.fastq.gz /scratch/gatoo/HSV1-mouse/"$name"_R.fastq.gz

      echo Starting mapping of "$name" against the virus
      echo Starting mapping sample "$name" against the virus with bowtie2

       bowtie2 -t -x ../HSV1-mouse/HSV1 -1 /scratch/gatoo/HSV1-mouse/"$name"_F.fastq.gz -2  /scratch/gatoo/HSV1-mouse/"$name"_R.fastq.gz -S mapped_"$name".sam
       samtools flagstat mapped_"$name".sam > flagstat_"$name".txt
       samtools view -h -F 4 mapped_"$name".sam | samtools view -bS > mapped_"$name".bam
       #Look at the alignment using samtools
       samtools view -b mapped_"$name".bam | samtools fillmd -e - /home/gatoo/HSV1-mouse/HSV1.fa | more

echo Finishing sample "$name"
      #tophat2 --num-threads 32 -o /scratch/gatoo/HSV1-mouse/output_HSV1_"$name" /scratch/gatoo/HSV1-mouse/HSV1 /scratch/gatoo/HSV1-mouse/"$name"_F.fastq.gz /scratch/gatoo/HSV1-mouse/"$name"_R.fastq.gz

    echo Finishing mapping "$name"
    
    echo Starting cufflinks analysis of "$name"
    
cufflinks --num-threads 32 -o /scratch/gatoo/HSV1-mouse/output_cuff_HSV1_"$name" --multi-read-correct -G /scratch/gatoo/HSV1-mouse/HSV1.gff3 /scratch/gatoo/HSV1-mouse/output_HSV1_"$name"/accepted_hits.bam
cufflinks --num-threads 32 -o /scratch/gatoo/HSV1-mouse/output_cuff_"$name" --frag-bias-correct /scratch/gatoo/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa --multi-read-correct -G /scratch/gatoo/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf /scratch/gatoo/HSV1-mouse/output_"$name"/accepted_hits.bam

cp  /scratch/gatoo/HSV1-mouse/output_cuff_"$name"/genes.fpkm_tracking "$name".txt

   echo $name counting: done!

done

echo Running cuffmerge on all samples

######################################################
#####Make sure you have an assemblies.txt file!! #####
######################################################

cuffmerge -g /scratch/gatoo/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf -s /scratch/gatoo/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.fa -p 8 assemblies.txt
#cuffmerge -g /scratch/gatoo/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf -s /scratch/gatoo/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.fa -p 8 assemblies-Male.txt

##Run Cuffdiff by using the merged transcriptome assembly along with the BAM files from TopHat for each replicate
#cuffdiff -o diff_out -b /scratch/gatoo/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.fa -p 8 -L "$name"-dLAT,"$name"-U1 -u merged_asm/merged.gtf \/scratch/gatoo/HSV1-mouse/output_"$name"-dLAT/accepted_hits.bam  \/scratch/gatoo/HSV1-mouse/output_"$name"-U1/accepted_hits.bam
#cuffdiff --num-threads 32 -o diff_out_"$name" -b /scratch/gatoo/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.fa -p 8 -L "$name"-dLAT,"$name"-U1 -u merged_asm/merged.gtf /scratch/gatoo/HSV1-mouse/output_"$name"-dLAT/accepted_hits.bam /scratch/gatoo/HSV1-mouse/output_"$name"-U1/accepted_hits.bam

###Run this to get the gene expression Male vs Female Specific condition (INDIVIDUALS) in this case dLAT vs WT
names='Fem Male'
for name in $names
 do
 echo Running cuffmerge on $name
  cuffdiff --num-threads 32 -o diff_out_2_"$name" -b /scratch/gatoo/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.fa -p 8 -L "$name"-dLAT,"$name"-WT -u merged_asm/merged.gtf /scratch/gatoo/HSV1-mouse/output_"$name"-dLAT/accepted_hits.bam /scratch/gatoo/HSV1-mouse/output_"$name"-WT/accepted_hits.bam
 done
 
cuffdiff --num-threads 32 -o diff_out_test -b /scratch/gatoo/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.fa -p 8 -L Fem-WT,Fem-U1,Fem-dLAT,Male-WT,Male-U1,Male-dLAT -u merged_asm/merged.gtf /scratch/gatoo/HSV1-mouse/output_Fem-WT/accepted_hits.bam /scratch/gatoo/HSV1-mouse/output_Fem-U1/accepted_hits.bam /scratch/gatoo/HSV1-mouse/output_Fem-dLAT/accepted_hits.bam /scratch/gatoo/HSV1-mouse/output_Male-WT/accepted_hits.bam /scratch/gatoo/HSV1-mouse/output_Male-U1/accepted_hits.bam /scratch/gatoo/HSV1-mouse/output_Male-dLAT/accepted_hits.bam
#done


echo All done
