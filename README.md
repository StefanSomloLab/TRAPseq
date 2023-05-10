# A TRAPseq protocol and analysis pipeline to identify mutant-cell specific transcriptome
Targeted purification of polysomal mRNA (TRAP-Seq) is a recently developed transcriptome profiling technique that maps translating mRNAs under specific tagged-cells. Herein, we describe our TRAPseq protocol for transcriptome profiling. Our protocol outlines all necessary steps for transcriptome profiling including the bioinformatic pipeline steps required to process, analyze, and identify differentially expressed genes (DEGs). <BR>
see, https://www.google.com/ for more TRAPseq experiment details.

<img src="https://github.com/macroant/TRAPseq/blob/main/doc/overview.png" height="400" width="700">
## 

To replicate results presented in the manuscript above, please run the script after installation of all required libraries.


##Install required libraries
We are using "R version 3.4.1" here.
For R, 'DESeq2' library is required, which can be installed via bioconductor https://bioconductor.org/packages/release/bioc/html/DESeq2.html
```R
> if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

> BiocManager::install("DESeq2")
```

## Workflow
The bioinformatics pipeline used for TRAPseq analysis consists of six main steps. Each of these steps performs a required function in processing the raw sequencing files (Fastq files) to obtain high-quality DEG list while the optional enrichment step performs the biological interpretation of DEGs once they are called. The six bioinformatics analysis steps are detailed below.<BR>
    
## <a name="workflow"></a> Step 1 Trimming & Alignment
Reads quality trimming and cleaned sequencing reads alignment to reference genomes. 
Firstly, [fastp](https://github.com/OpenGene/fastp) is run to remove low quality bases from raw sequencing reads while also removing any potential adapter sequences. <BR>
```bash
$ fastp -i R1.fastq.gz -I R2.fastq.gz \
    -o R1.fastq.gz.out.fq.gz -O R2.fastq.gz.out.fq.gz \
    -j test.json -h test.html -w threads_num
```
fastp can be generated detailed HTML reports to assess reads stats, quality, GC content, duplication rate, length distribution, K-mer content and adapter contamination.<BR>
Once the reads are cleaned, alignment to the reference genomes is performed using [STAR](https://github.com/alexdobin/STAR). In this study, reads are aligned to the mouse genome by STAR 2.7.9. Following alignment, unaligned reads are removed, and alignment statistics are calculated.<BR>
#prior to alignment, STAR require you to construct and index the genome.
```bash
$ STAR --runThreadN threads_num --runMode genomeGenerate --genomeDir /path-to-index \
		 --genomeFastaFiles /path-to-fasta \
		 --sjdbGTFfile /path-to-GTF-file \
		 --sjdbOverhang 99
```
#alignment with STAR
```bash
$ STAR --genomeDir /path-to-STAR-index --readFilesIn R1.fastq.gz.out.fq.gz R2.fastq.gz.out.fq.gz \
    --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 16000000000 --outSAMunmapped Within \
    --runThreadN threads_num --outFileNamePrefix output_test
```

## <a name="workflow"></a> Step 2 Raw Reads Quantification & Normalization
Assign aligned reads to genomic features (genes) and summarize raw reads count for each genomic feature and each sample. Normalize raw read counts to remove technical biases, e.g. sequencing depth and gene length, and make normalized gene expression values directly comparable within and across samples.
The gene abundance quantification tool [featureCount](https://subread.sourceforge.net/featureCounts.html) is used to calculate raw read count for each gene and each sample. 

```bash
$ featureCounts -T threads_num -p -t exon -g gene_id \
		-a  /path-to-GTF-file \
		-o featureCounts.gene_level_exon_counts_PE_samples.txt \
		-J *.bam
```
We used the gene expression detection threshold where genes were selected if there were ≥6 reads in at least 20% of samples. <BR>
```r
#r function for filtering out low/non-expressed genes
#x is the raw reads counts for a gene across samples
filter_function <- function (x) {
  y <- as.numeric(x)
  if ((sum(y>=6)/length(y))>=0.2) {
    return(TRUE)
  }else{
    return(FALSE)
  }
}
```
The relative abundance of a transcript is calculated by Fragments per Kilobase per Million Mapped Fragments (FPKM) and Transcripts per million (TPM).<BR>
```r
#r function for TPM normalization
#x is the raw reads counts for genes in a sample, y is their corresponding length of genes

tpm_Normalization <- function(x, y) {
  rate <- x / y
  (rate / sum(rate) ) * 1e6
}

```
## <a name="workflow"></a> Step 3 Detection of Outlier Samples
Evaluate and remove outlier samples can significantly improve the reliability of DEGs and downstream functional analysis. We used [robustPCA](https://cran.r-project.org/web/packages/rrcov/index.html) for outlier detection as it was found to result in stable performance and  the default cutoff value allows correct identification of outliers without yielding false positives in an evaluation study.<BR>

```r
# data_matrix_log, columns are samples, rows are genes, expression values are log2 transformed TPM values
> pcaData <- prcomp(t(data_matrix_log))
> pcaVar <- (pcaData$sdev)^2/sum((pcaData$sdev)^2) * 100

# Scree Plot
> plot(pcaVar, xlab = "Component Number", 
     ylab = "% of variance", 
     main = "scree plot")
> abline(v=3, lty = 2, col = "red")
> text(x=5, y=20, paste("top3 = ", round(sum(pcaVar[1:3]),1), "%", sep=""), col = "red")
```
<img src="https://github.com/macroant/TRAPseq/blob/main/doc/scree_plot.png" height="400" width="500">
	
```r
# Robust PCA, k=3
> pc <- PcaGrid(t(data_matrix_log), k = 3, method = "qn")
> plot(pc, main = "Robust PCA, k=3", xlim = c(0, 3.5))
```
<img src="https://github.com/macroant/TRAPseq/blob/main/doc/robustPCA.png" height="400" width="500">
<BR>
#Visualization TRAPseq Samples with classic PCA
<img src="https://github.com/macroant/TRAPseq/blob/main/doc/PCA.png" height="400" width="500">
<BR>

## <a name="workflow"></a> Step 4 Cell Type Deconvolution
TRAPseq enriches specific cells than other cells, cell type deconvolution analysis can help reveal the cell type composition of TRAPseq transcriptome profiles. We estimated the cellular compositions of the TRAPseq samples through deconvolution analysis. We obtained gene expression signatures for major tubules segments from a mouse micro-dissected tubules dataset (GSE150338, PMID: 33769951), and used as input into [CIBERSORTx](https://cibersortx.stanford.edu) to estimate the cellular composition of our TRAPseq data. <BR>

## <a name="workflow"></a> Step 5 Detection of Differential Expressed Genes & Filtering
This step detects genes with statistically significant differences or changes in read counts or expression levels between two experimental conditions. The statistical significance test is conducted using the [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) with negative binomial distribution and Benjamini-Hochberg adjusted p-values (FDR). To screen for more meaningful DEGs, we only consider genes that are protein-coding genes and have the same expression change direction in comparison of SKO/WT and SKO/DKO. <BR>

```r
#Detection of DEGs
# matrix is raw reads count matrix, columns are samples, rows are genes
##                       Sample1 Sample2 Sample3 ...
## ENSMUSG00000033845.13     4565     4529     3560 ...
## ENSMUSG00000025903.14     7539     7552     5545 ...
	...                   ...      ...      ...
# phe is sample information table
##    SampleID Pheno Sex
## 2   Sample1    WT   M
## 3   Sample2    WT   M
## 9   Sample3    WT   M
## 10  Sample4   SKO   M
## 11  Sample5   SKO   M

> DES_results <- DESeqDataSetFromMatrix(countData = matrix,
                                      colData = phe,
                                      design = ~ Pheno)
> DES_results <- relevel(DES_results$Pheno, ref="WT")
> DES_results <- DESeq(DES_results)
> DES_results <- results(DES_results)
> DES_results <- as.data.frame(DES_results)
> DES_results$GeneID <- row.names(DES_results)
```
#r code for making valcano plot			
```r
# SKO_WT, SKO_DKO, and DKO_WT are 3 DEGs lists generated by DESeq2
SKO_WT <- na.omit(SKO_WT)
SKO_WT$SigCode <- "NoSig"
#SKO_WT[SKO_WT$padj<=0.25,]$SigCode = "FDR0.25"
nrow(SKO_WT[SKO_WT$padj<=0.25,])
nrow(SKO_WT[SKO_WT$padj<=0.05,])
SKO_WT[SKO_WT$padj<=0.05 & SKO_WT$log2FoldChange < 0,]$SigCode = "Down"
SKO_WT[SKO_WT$padj<=0.05 & SKO_WT$log2FoldChange >= 0,]$SigCode = "Up"
SKO_WT$SigCode <- factor(SKO_WT$SigCode, levels = c("Up", "Down", "NoSig"), ordered = T)

SKO_DKO <- na.omit(SKO_DKO)
SKO_DKO$SigCode <- "NoSig"
#SKO_DKO[SKO_DKO$padj<=0.25,]$SigCode = "FDR0.25"
nrow(SKO_DKO[SKO_DKO$padj<=0.25,])
nrow(SKO_DKO[SKO_DKO$padj<=0.05,])
SKO_DKO[SKO_DKO$padj<=0.05 & SKO_DKO$log2FoldChange < 0,]$SigCode = "Down"
SKO_DKO[SKO_DKO$padj<=0.05 & SKO_DKO$log2FoldChange >= 0,]$SigCode = "Up"
SKO_DKO$SigCode <- factor(SKO_DKO$SigCode, levels = c("Up", "Down", "NoSig"), ordered = T)

DKO_WT <- na.omit(DKO_WT)
DKO_WT$SigCode <- "NoSig"
#DKO_WT[DKO_WT$padj<=0.25,]$SigCode = "FDR0.25"
nrow(DKO_WT[DKO_WT$padj<=0.25,])
nrow(DKO_WT[DKO_WT$padj<=0.05,])
DKO_WT[DKO_WT$padj<=0.05 & DKO_WT$log2FoldChange < 0,]$SigCode = "Down"
DKO_WT[DKO_WT$padj<=0.05 & DKO_WT$log2FoldChange > 0,]$SigCode = "Up"
DKO_WT$SigCode <- factor(DKO_WT$SigCode, levels = c("Up", "Down", "NoSig"), ordered = T)

#making volcano plot
SKO_WT_plot <- ggplot(data = SKO_WT, aes(x = log2FoldChange, y = -log10(padj), col = SigCode, label = GeneName) ) + 
  scale_color_manual(values = c("red", "blue", "gray")) + theme_light() + 
  geom_point() + theme(legend.position = "none", 
           panel.border = element_blank(), 
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(), 
           axis.line = element_line(colour = "black")) + 
  geom_point(data = subset(SKO_WT, show == "Y"), aes(x = log2FoldChange, y = -log10(padj)), color="black") + 
  geom_text_repel(data = subset(SKO_WT, show == "Y"), aes(label = GeneName), 
                  color="black", size = 5) +
  geom_vline(xintercept = c(log2(0.5), log2(2)), col = "gray", linetype = 2) + 
  geom_hline(yintercept = c(-log10(0.05)), col = "red", linetype = 2) + 
  labs(title = "SKO vs WT, M") + xlim(-5, 5) + scale_y_sqrt(breaks = c(0, 30, 60, 90, 120), limits = c(0, 120)) 
  

SKO_DKO_plot <- ggplot(data = SKO_DKO, aes(x = log2FoldChange, y = -log10(padj), col = SigCode, label = GeneName) ) + 
  scale_color_manual(values = c("red", "blue", "gray")) + theme_light() + 
  geom_point() + theme(legend.position = "none", 
           panel.border = element_blank(), 
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(), 
           axis.line = element_line(colour = "black")) + 
  geom_point(data = subset(SKO_DKO, show == "Y"), aes(x = log2FoldChange, y = -log10(padj)), color="black") + 
  geom_text_repel(data = subset(SKO_DKO, show == "Y"), aes(label = GeneName), 
                  color="black", size = 5) +
  geom_vline(xintercept = c(log2(0.5), log2(2)), col = "gray", linetype = 2) + 
  geom_hline(yintercept = c(-log10(0.05)), col = "red", linetype = 2) + 
  labs(title = "SKO vs DKO, M") + xlim(-5, 5) + scale_y_sqrt(breaks = c(0, 30, 60, 90, 120), limits = c(0, 120))

DKO_WT_plot <- ggplot(data = DKO_WT, aes(x = log2FoldChange, y = -log10(padj), col = SigCode, label = GeneName) ) + 
  scale_color_manual(values = c("red", "blue", "gray")) + theme_light() + 
  geom_point() + theme(panel.border = element_blank(), 
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(), 
           axis.line = element_line(colour = "black")) + 
  geom_point(data = subset(DKO_WT, show == "Y"), aes(x = log2FoldChange, y = -log10(padj)), color="black") + 
  geom_text_repel(data = subset(DKO_WT, show == "Y"), aes(label = GeneName), 
                  color="black", size = 5) +
  geom_vline(xintercept = c(log2(0.5), log2(2)), col = "gray", linetype = 2) + 
  geom_hline(yintercept = c(-log10(0.05)), col = "red", linetype = 2) + 
  labs(title = "DKO vs WT, M") + xlim(-5, 5) + scale_y_sqrt(breaks = c(0, 30, 60, 90, 120), limits = c(0, 120)) 

SKO_WT_plot + SKO_DKO_plot

combined_plot <- SKO_WT_plot + SKO_DKO_plot + DKO_WT_plot

combined_plot
```
<img src="https://github.com/macroant/TRAPseq/blob/main/doc/volcano.png" height="400" width="700">

## <a name="workflow"></a> Step 6 Enrichment Analysis (optional analysis step)
Biological interpretations of selected DEGs are performed using various bioinformatics tools, including [DAVID](https://david.ncifcrf.gov), [Metascape](https://metascape.org/gp/index.html#/main/step1), and commercal tool [IPA](https://digitalinsights.qiagen.com/products-overview/discovery-insights-portfolio/analysis-and-visualization/qiagen-ipa/) or [MetaCore](https://clarivate.com/products/biopharma/discovery-clinical-regulatory/early-research-intelligence-solutions/).

## Useful References
