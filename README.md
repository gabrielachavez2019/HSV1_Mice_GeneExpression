# GeneExpression Analysis of the HSV1 infected mice 2021
This is a repository for gene expression analysis of the RNA-Seq data got by Ilumina from infected mice with HSV1

## RNA-Seq PIPELINE
This pipeline performs the following tasks:
- perform quality control on FastQ files (using FastQC)
- align reads of each sample against reference genome (using tophat2) 
- count reads in features (cufflinks)
- normalize read counts (using cufflinks::cuffdiff)
- calculate RPKMs (using CummeRbund)
- perform DE analysis for standard designs (using CummeRbund)
- Generate Gene Onthology and WikiPathways outputs!

## Download
Use git clone:
```
git clone git@github.com:gabrielachavez2019/HSV1_Mice_GeneExpression.git
```

## Installation
Download the RNAseq pipeline.


## Genome files
#### Generate genome indexes or for Mice, just download it from [here](https://uswest.ensembl.org/Mus_musculus/Info/Index)
Create a directory where you want to store the indexes (e.g. /GENOMES/Mus_musculus/UCSC/mm10/).

Collect the following files for your genome:
- Fasta file containing the genome reference sequences
- GTF file containing annotated transcripts

#### Generate genome indexes
Generate genome indexes files for the Viral genomes. 
In this case we are using the Human herpesvirus 1 strain Mckrae, available at the GenBank:[JX142173.1](https://www.ncbi.nlm.nih.gov/nuccore/JX142173). 
Get the FASTA file from [HSV1.FASTA](https://www.ncbi.nlm.nih.gov/nuccore/JX142173.1?report=fasta)
The genome indexes are saved to disk and need only be generated once for each genome/annotation combination.

```  
bowtie2-build HSV1.FASTA HSV1
```

## Usage
#### Run pipeline
```bash
sbatch mapping.sh -input [/path/to/rundir] -outputDir [/path/to/outputdirname] -mail [email]
```

## Dependencies
#### Load modules
- fastqc 0.11.7 
- bowtie2 2.3.4.1  
- samtools 1.10
- trimmomatic 0.38
- cufflinks 2.2.1
- tophat 2.1.2 
- R 4.0.3    

#### R packages
- DESeq2 1.6.3
- edgeR 3.8.6
- ggplot2
- gplots
- RColorBrewer
- CummeRbund


#### Mapping
Illumina pair-end hight quality reads are aligned first to the HSV1 genome using minimap, which is better for viral genomes. Followed by mapping them agaist the mice genome using TOPHAT2, which was designed to specifically address many of the challenges of RNA-seq data mapping, and uses a novel strategy for spliced alignments.

#### Read counting
Counting sequencing reads in features (genes/transcripts/exons) is done with cufflinks using the union mode.

#### Normalizing read counts
Cufflinks with default options --library-norm-method classic-fpkm (default) for Cufflinks --library-norm-method geometric (default) for cuffdiff.
Cuffnorm will report both FPKM values and normalized, estimates for the number of fragments that originate from each gene, transcript, TSS group, and CDS group. Note that because these counts are already normalized to account for differences in library size, they should not be used with downstream differential expression tools that require raw counts as input.

#### Calculate FPKMs
Cuffdiff calculates the FPKM of each transcript, primary transcript, and gene in each sample. Primary transcript and gene FPKMs are computed by summing the FPKMs of transcripts in each primary transcript group or gene group. The results are output in FPKM tracking files in the format described here. There are four FPKM tracking files:

isoforms.fpkm_tracking	Transcript FPKMs
genes.fpkm_tracking	Gene FPKMs. Tracks the summed FPKM of transcripts sharing each gene_id
cds.fpkm_tracking	Coding sequence FPKMs. Tracks the summed FPKM of transcripts sharing each p_id, independent of tss_id
tss_groups.fpkm_tracking	Primary transcript FPKMs. Tracks the summed FPKM of transcripts sharing each tss_id

#### Differential expression analysis
CummeRbund makes managing and querying data easier by loading the data into multiple objects of several different classes and having functions to query them. Because all of this gets stored in an sql database, you can access it quickly without loading everything in to memory.

  readCufflinks- Most important function designed to read all the output files that cuffdiff generates into an R data object (of class CuffSet).
 ```
 cuff_data <- readCufflinks(diff_out)
``` 
CummeRbund has at least 6 classes that it will place different parts of your data in.
Now you can access information using different functions:  gene information using genes(cuff_data), your isoform level output using isoforms(cuff_data), TSS related groups using tss(cuff_data) and so forth

You can explore global statistics on data for quality, dispersion, distribution of gene expression scores etc or Plot things for specific features or genes.

#### Additional tools
Please contact Gabriela Toomer (gabriela.toomer@okstate.edu) if you want to add additional tools/scripts/options or have any questions.
