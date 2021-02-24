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
```bash
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

`` bash 
bowtie2-build HSV1.FASTA HSV1
``

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
RNA-seq reads are aligned to the reference genome using STAR, which was designed to specifically address many of the challenges of RNA-seq data mapping, and uses a novel strategy for spliced alignments.

#### Detecting fusion genes
The pipeline searches for fusion genes / chimeric junctions by default. The default minimum segment length is 1 and can be edited with option -chimSegmentMin.

#### Read counting
Counting sequencing reads in features (genes/transcripts/exons) is done with htseq-count using the union mode.
The raw read counts table can be found in the directory: <rundir>/read_counts/<run>_raw_counts.txt. This table contains the Ensembl (gene/transcript/exon) IDs (rows) and the raw read 
counts per library (columns).

#### Normalizing read counts
The raw read counts are normalized using the DESeq method included in the DESeq Bioconductor package and is based on the hypothesis that most genes are not DE. A DESeq scaling factor for
a given lane is computed as the median of the ratio, for each gene, of its read count over its geometric mean across all lanes. The underlying idea is that non-DE
genes should have similar read counts across samples, leading to a ratio of 1. Assuming most genes are not DE, the median of this ratio for the lane provides an
estimate of the correction factor that should be applied to all read counts of this lane to fulfill the hypothesis. By calling the estimateSizeFactors() and sizeFactors()
functions in the DESeq Bioconductor package, this factor is computed for each lane, and raw read counts are divided by the factor associated with their sequencing lane. 

#### Calculate RPKMs
RPKM is a method of quantifying gene expression from RNA sequencing data by normalizing for total read length and the number of sequencing reads. RPKMs are calculated using the RPKM()
function included in the edgeR Bioconductor package.
RPKM = [# of mapped reads]/([length of transcript]/1000)/([total reads]/10^6)

CAUTION: Make sure the RPKM normalization is actually necessary in your analysis. If you are comparing gene expression among samples only, there really is no reason to normalize by length
as you will be dividing each gene among the samples by a constant (gene length). You only need to use RPKM when you are comparing transcript expression within one sample.

#### Differential expression analysis
Differential expression analysis can be done only for standard designs, such as:
   sample1	test
   sample2	control
   sample3	test
   sample4	test
   sample5	control

DE analysis is done using the DESeq2 Bioconductor package. It takes the merged raw read counts (from HTseq-count) as an input and creates the following output inside the /<run>/DEanalysis folder:
- <run>_DEanalysis_all.txt : result table of DE analysis (ordered by p-value). The columns are:
    - Ensembl (gene/transcript/exon) ID
    - baseMean -> the average of the normalized count values, dividing by size factors, taken over all samples
    - log2FoldChange -> effect size estimate (tells us how much the gene's expression seems to have changed due to for example treatment in comparison to untreated samples
    - lfcSE -> standard error estimate for the log2 fold change estimate
    - stat
    - pvalue -> indicates the probability that a fold change as strong as the observed one, or oven stronger, would be seen under the situation described by the null hypothesis.
    - padj -> adjusted p-value using the Benjamini-Hochberg adjustment
    - gene_id
    - gene_name
- <run>_MAplot.png: shows the log2 fold changes attributable to a given variable over the mean of normalized counts. Points are colored red if the adjusted p value is less than 0.1, points that fall out of the window are plotted as open triangles pointing either up or down.
- <run>_sampletosample_distances.png: heatmap of the sample-to-sample distances.
- <run>_PCAplot.png: principal component plot of the samples.


#### Additional tools
Please contact Gabriela Toomer (gabriela.toomer@okstate.edu) if you want to add additional tools/scripts/options or have any questions.
