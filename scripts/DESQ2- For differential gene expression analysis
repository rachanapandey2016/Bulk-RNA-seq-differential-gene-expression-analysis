library(DESeq2)
#tximport is an R/Bioconductor package designed to convert transcript-level quantification 
#data (e.g. from Salmon, Kallisto, or RSEM) into gene-level counts or abundances for downstream #differential #expression analysis (e.g., with DESeq2, edgeR)
#Many RNA-seq tools (e.g., DESeq2) work best with gene-level counts. tximport aggregates 
#transcript-level data into gene-level summaries. it Merge multiple transcripts of same gene into 
#single expression value.
library(tximport)  #input file is quant.sf which has abundance estimates
library(data.table)
library(ashr)
library(ggplot2) #for plotting
library(ggrepel)
setwd("C:/PhD_Courses/Fourth Semester/GCD 8141 Computational genomics/Project-3/Count matrix/S4/")
# Read in the gene/transcript IDs
#The tximport function needs two file .sh file from the salmon and annotation file which converts #transcript id to gene id. tx2gnene is transcript-to-gene mapping table, transcript ID is Unique #identifier for each transcript and ensgene is Unique identifier for the gene. it is used by tximport to #group transcripts from the same gene.
tx2gene = fread("./tx2gene_grch38_ens94.txt",sep = "\t")
head(tx2gene)

#I already have all 12 samples quantification files (quant.sf from Salmon) in 
#a single folder, and now I will compile them into a single count matrix.
#First List All quant.sf Files
files <- list.files(pattern = "quant.sf", full.names = TRUE, recursive = TRUE)
files

#Extract sample names (remove "./" and ".quant.sf")
sample_names <- gsub("^\\./(.*)\\.quant\\.sf$", "\\1", files)

#Assign names
names(files) <- sample_names
names(files)

#With both files, let’s load the transcript counts and convert them to gene counts.
txi <- tximport(files, 
                type = "salmon", 
                tx2gene = tx2gene [, c("tx_id","ensgene")], 
                countsFromAbundance = "lengthScaledTPM",
                ignoreTxVersion = TRUE)  # Set FALSE if transcript IDs contain versions

#Salmon outputs gene expression in TPM (transcript per million) units, which normalizes for read depth and gene length. DESeq2 requires its own normalization method 
#and therefore, requires un-normalized gene counts. tximport will automatically do this when we read in the TPM values.  
#Now, to check if your genes are differentially expressed between two treatments, we need to provide the treatment level of each sample. 

#I have already extracted the sample names earlier and store it in sample_names variable.
#Now we can directly create metadata
meta <- data.frame(
  sample_id = sample_names,  # Use the already-cleaned names
  group = factor(gsub("_.*", "", sample_names)),  # Extract prefix (e.g., "S4L")
  row.names = sample_names  # Set sample names as row IDs
)
meta

#Checking that sample names match in both files
all(colnames(txi$counts) %in% rownames(meta))   #Should return TRUE

#To use DESeq for DGE, we need to store our data in DESeqDataSet format. This is a convenient data structure made up of lists and data.frames and we will look at some 
#of the things stored in it. It requires three arguments: (i) the un-normalized gene counts loaded with tximport, (ii) a data.frame that contains information about the 
#treatment for each sample (or the group the sample belongs in, in the same order as the columns in the count data), and (iii) a design formula that tells DESeq how to 
#model the predictors on gene counts.

#Creating DESeq2 object
dds <- DESeqDataSetFromTximport(
  txi = txi,
  colData = meta,
  design = ~ group  # Using the 'group' column as it has (S4L vs S4S)
)

#Now we will do PCA to look at batch effects.
# Regularized log transformation (stabilizes variance)
rld <- rlog(dds, blind = TRUE)

pca_plot<-plotPCA(rld, intgroup = "group")

# Save the PCA plot as a PDF file
ggsave("pca_plot.pdf", plot = pca_plot, width = 6, height = 6)

pca_plot

# 1. Get PCA results from rlog-transformed data
rld_mat <- assay(rld)  # Extract rlog counts
pca <- prcomp(t(rld_mat))  # Transpose for samples-as-rows

# 2. Combine PCA results (PC3 & PC4) with metadata
df <- data.frame(meta, pca$x[,3:4])  # Select PC3 and PC4
colnames(df)[(ncol(meta)+1):(ncol(meta)+2)] <- c("PC3", "PC4")  # Rename columns
ggplot(df, aes(x = PC3, y = PC4, color = group)) +  # Use 'group' from your meta
  geom_point(size = 3) +
  labs(title = "PC3 vs PC4 -Stage 4 Melanoma (short progressor vs Long responder",
       x = paste0("PC3 (", round(100*pca$sdev[3]^2/sum(pca$sdev^2), 1), "% variance)"),
       y = paste0("PC4 (", round(100*pca$sdev[4]^2/sum(pca$sdev^2), 1), "% variance)")) +
  theme_minimal() +
  scale_color_manual(values = c("S4L" = "blue", "S4S" = "red")) 


#Run analysis
dds <- DESeq(dds)

#Defining contrast
contrast_melanoma <- c("group", "S4L", "S4S")  # Short Progressors (S4L) vs S4S)

# Extract results (with FDR cutoff of 0.05)
res <- results(dds, contrast = contrast_melanoma, alpha = 0.05)

#view the top differentially expressed genes
head(res)

plotMA(res, ylim = c(-2, 2), colSig = "red", main = "MA Plot without shrinkage in long responder vs short progressor")

#Applyfoldchangeshrinkage
 res_tableOE_shrunk<-lfcShrink(dds,contrast=contrast_melanoma,
 res=res,type="ashr")

 plotMA(res_tableOE_shrunk,ylim=c(-2,2),colSig="red", main = "MA Plot after shrinkage in long responder vs short progressor")


 ##foldchanges
contrast_mel_shr<-c("group", "S4L", "S4S")
res_tableKD<-results(dds,contrast=contrast_mel_shr,alpha=0.05)
res_shrunk<-lfcShrink(dds,contrast=contrast_mel_shr,res=res_tableKD,type="ashr")


summary(res_tableOE_shrunk)

# Add gene IDs as a column (these are currently the rownames)
res_df <- as.data.frame(res_shrunk)
res_df$gene <- rownames(res_df)

# Select only the desired columns
res_subset <- res_df[, c("baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj", "gene")]

# Save to CSV
#write.csv(res_subset, file = "Stage4_S4L_vs_S4S_DEGs_results.csv", row.names = FALSE)


#Removing the duplicated gene
tx2gene_unique <- tx2gene[!duplicated(tx2gene$ensgene), ]

res_annotated <- merge(res_df, tx2gene_unique, by.x = "gene", by.y = "ensgene")

# calculate significance column
res_annotated$sig <- ifelse(!is.na(res_annotated$padj) & res_annotated$padj < 0.05, "s", "ns")
 
# Filter to all significant genes with valid gene symbols
sig_genes <- subset(res_annotated, padj < 0.05 & !is.na(symbol) & symbol != "")

# Check how many are being labeled
nrow(sig_genes)  # Should print 15
S4l_DEGs<- ggplot(res_annotated) +
  geom_point(aes(log2FoldChange, -log10(pvalue), color = sig), alpha = 0.6) +
  geom_text_repel(data = sig_genes,
                  aes(log2FoldChange, -log10(pvalue), label = symbol),
                  size = 3,
                  max.overlaps = Inf,    # ensures ALL labels can be drawn
                  box.padding = 0.5,
                  segment.color = "grey50") +
  theme_classic() +
   theme( axis.title = element_text(face = "bold", size = 11),
          plot.title = element_text(face = "bold", size = 12)) +
  labs(title = "Volcano Plot: Differentially expressed gene in Stage4 Melanoma Long responder \n vs Short Progressor after ICI  Immunotherapy",
       x = "Log2 Fold Change",
       y = "-Log10(p-value)",
       color = "Significance")
ggsave("Volcano_Plot_for_S4l_DEGs.pdf", plot = S4l_DEGs, width = 8, height = 6)
S4l_DEGs

