# DEDSeq
A pipeline to study **d**ifferential gene **e**xpression and **d**ifferential **s**plicing using RNA**seq** data

# Introduction: 
DEDSeq is a pipeline to study differential gene expression and differential splicing using RNAseq data. 

Operation system: Linux  
Languages used: bash, Perl, R  
R packages needed: DESeq2, ggplot2, edgeR, BiocParallel, fmsb, gplots, openxlsx, DEXSeq, knitr, pander, rmarkdown  

Other tools needed:   
STAR  
samtools  
bedtools  

# Main functions:
1, map the RNAseq reads in fastq format to a genome using STAR, creating .bam files. Generate mapping index (.bai) files. Generate bigwig (.bw) files for visualization of read mapping profile.

2, quantify the gene expression based on the read counts mapped (in default, only CDS region of a protein-coding gene is used for read counting to reduce the effect of alternative polyadenylation in the quantification of gene expression);

3, study differential expression (DE) using DESeq2 (preferred method) or Fisher’s exact test (when no replicates are available);

4, generate an Excel table and a tab-delimited file of genes with their read counts, RPKM/FPKM values, log2 fold changes and P-values (from DEseq2) information;

5, generate a html format RNAseq-report for read mapping statistics, DE analysis. Generate plots: dentroplot, volcano plot, MA plot, scatter plot etc;

6, identify novel junctions from the read mapping, separate each gene to blocks (exons or intron regions) based on the junction reads and quantify read counts for each gene block;

7, study alternative splicing using gene-block expression using DEXSeq or Fisher’s exact test

8, study splicing using percent-spliced-in (PSI) method, can use t-test, Wilcoxon rank-sum test and/or Fisher’s exact test. t-test and Wilcoxon rank-sum test can only be used when there are enough replicates per sample group (recommend at least 3 for each group). The output files include tab-delimited files and an Excel table with information of read counts for inclusion/exclusion, PSI values, deltaPSI values, P-values, and regulation types.


# How to use:
1, copy the whole package to a linux server;

2, download or create necessary files in PIpeline/ReferenceDB folder. eg. the genome sequence files, STAR mapping index files etc. The final gene annotation files for human hg19 and mm9 are already in the packages. For other species or genome version, go to each sub-folder and download or create the annotation files using the bash commands, perl or R script provided in the sub-folder;

3, create a project folder for a specific RNAseq dataset, copy the “pipeline_template.sh” file into it and edit it based on the settings of the study. Eg. specify the genome version (eg. hg19, or mm9), whether the RNAseq is paired-end or single-end, stranded or not, sample names, sample group names and replicates information, and test and control samples etc. Follow the instruction in “pipeline_template.sh” file;

4, make pipeline_template.sh executable and run ./pipeline_template.sh in the shell or submit the job using a job scheduler in the cluster;

5, not all analysis is required. eg, if you already have mapped files in bam format, you can skip the mapping step.


