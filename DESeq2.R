#DESeq2 Analysis

# Load libraries
library("DESeq2")
library("apeglm")
library("pheatmap")
library("tidyverse")
library("ggrepel")
library("ashr")
library("goseq")
library("biomaRt")


# create an object with the directory containing your counts:

directory <- "../../vera_licona_counts/"

# ensure the count files are where you think they are
list.files(directory)

samples <- list.files(directory, pattern = ".*counts$")
meta <- read.csv("../../MetaData.csv")
# ensure that sampleFiles and metadata table are in the same order


all( str_remove(samples, ".counts") == meta[,1] )

# now create a data frame with sample names, file names and treatment information. 
sampleTable <- data.frame(
  sampleName = meta$Run,
  fileName = samples,
  age = meta$age_group
)
sampleTable
# create the DESeq data object
ddsHTSeq <- DESeqDataSetFromHTSeqCount(
  sampleTable = sampleTable, 
  directory = directory, 
  design = ~ age
)

######################################################
# Reset treatment factors
######################################################

# factor order determines reference level (i.e. control vs exposed)
# factors are automatically ordered alphabetically
# we want KC as the reference level

# To see the levels as they are now:
ddsHTSeq$age

# To replace the order with one of your choosing, create a vector with the order you want:
age <- c("young","old")

# Then reset the factor levels:
ddsHTSeq$age <- factor(ddsHTSeq$age, levels = age)

# verify the order
ddsHTSeq$age


###Statistical analyses

#expression across genes

#sumcounts for each gene
sumcounts <- rowSums(counts(ddsHTSeq))
#take the log
logsumcounts <- log(sumcounts,base=10)
hist(logsumcounts,breaks=100)


#Get results:

#genes that are differentially expressed in the pancreas in the different age groups (with young acting as control)
deseq <- DESeq(ddsHTSeq)
resultsNames(deseq)
results <- results(deseq, name="age_old_vs_young")
summary(results)


######################################################
# Get a table of shrunken log2 fold changes
######################################################

# get shrunken log fold changes
res_shrink <- lfcShrink(deseq,type="ashr",coef="age_old_vs_young")

# plot the shrunken log2 fold changes against the raw changes:

data.frame(l2fc=results$log2FoldChange, l2fc_shrink=res_shrink$log2FoldChange, padj=results$padj) %>%
  filter(l2fc > -5 & l2fc < 5 & l2fc_shrink > -5 & l2fc_shrink < 5) %>%
  ggplot(aes(x=l2fc, y=l2fc_shrink,color=padj > 0.1)) +
  geom_point(size=.25) + 
  geom_abline(intercept=0,slope=1, color="gray")


# get the top 20 genes by shrunken log2 fold change
# this will include genes with outliers
arrange(data.frame(res_shrink), -abs(log2FoldChange)) %>% head(., n=20)


######################################################
# Data visualization
######################################################

# MA plot
plotMA(results, ylim=c(-4,4))
plotMA(res_shrink, ylim=c(-4,4))

##############

#Volcano plot

# negative log-scaled adjusted p-values
log_padj <- -log(res_shrink$padj,10)
log_padj[log_padj > 100] <- 100

# plot
plot(x=res_shrink$log2FoldChange,
     y=log_padj,
     pch=20,
     cex=.5,
     col=(log_padj > 10)+1, # color padj < 0.1 red
     ylab="negative log-scaled adjusted p-value",
     xlab="shrunken log2 fold changes")

#############

# PCA plot

# normalized, variance-stabilized transformed counts for visualization
vsd <- vst(deseq, blind=FALSE)

plotPCA(vsd, intgroup=c("age"))

# alternatively, using ggplot

data <- plotPCA(vsd,returnData=TRUE,intgroup=c("age"))
data
p <- ggplot(data,aes(x=PC1,y=PC2,col=paste(age)))
p <- p + geom_point() + 
  xlab(paste("PC1: ", round(attr(data,"percentVar")[1],2)*100, "% variation explained", sep="")) + 
  ylab(paste("PC2: ", round(attr(data,"percentVar")[2],2)*100, "% variation explained", sep="")) +
  geom_label_repel(aes(label=name))
p

###############################################################################
########################## NEW ANALYSIS #######################################

####Remove outlier######

samples <- samples[!grepl("SRR1299340.counts",samples)]
samples
meta <- read.csv("../../MetaData.csv")
meta <- filter(meta,!grepl("SRR1299340", Run))
meta
all( str_remove(samples, ".counts") == meta[,1] )

# now create a data frame with sample names, file names and treatment information. 
sampleTable <- data.frame(
  sampleName = meta$Run,
  fileName = samples,
  age = meta$age_group
)
sampleTable

# create the DESeq data object
ddsHTSeq <- DESeqDataSetFromHTSeqCount(
  sampleTable = sampleTable, 
  directory = directory, 
  design = ~ age
)

######################################################
# Reset treatment factors
######################################################

# factor order determines reference level (i.e. control vs treatment)
# factors are automatically ordered alphabetically
# we want young as the reference level

# To see the levels as they are now:
ddsHTSeq$age

# To replace the order with one of your choosing, create a vector with the order you want:
age <- c("young","old")

# Then reset the factor levels:
ddsHTSeq$age <- factor(ddsHTSeq$age, levels = age)

# verify the order
ddsHTSeq$age


###Statistical analyses

#expression across genes

#sumcounts for each gene
sumcounts <- rowSums(counts(ddsHTSeq) > 10)

#take the log
logsumcounts <- log(sumcounts,base=10)
hist(logsumcounts,breaks=100)
#keep genes where 4 or more samples have a count more than 10
keep <- sumcounts >= 4
ddsHTSeq <- ddsHTSeq[keep,]
dds <- DESeq(ddsHTSeq)
resultsNames(dds)
results <- results(dds, name="age_old_vs_young")
summary(results)


######################################################
# Get a table of shrunken log2 fold changes
######################################################

# get shrunken log fold changes
res_shrink <- lfcShrink(dds,type="ashr",coef="age_old_vs_young")

# plot the shrunken log2 fold changes against the raw changes:

data.frame(l2fc=results$log2FoldChange, l2fc_shrink=res_shrink$log2FoldChange, padj=results$padj) %>%
  filter(l2fc > -5 & l2fc < 5 & l2fc_shrink > -5 & l2fc_shrink < 5) %>%
  ggplot(aes(x=l2fc, y=l2fc_shrink,color=padj > 0.1)) +
  geom_point(size=.25) + 
  geom_abline(intercept=0,slope=1, color="gray")


# get the top 20 genes by shrunken log2 fold change
# this will include genes with outliers
arrange(data.frame(res_shrink), -abs(log2FoldChange)) %>% head(., n=20)


######################################################
# Data visualization
######################################################

# MA plot
plotMA(results, ylim=c(-4,4))
plotMA(res_shrink, ylim=c(-4,4))

##############

#Volcano plot

# negative log-scaled adjusted p-values
log_padj <- -log(res_shrink$padj,10)
log_padj[log_padj > 100] <- 100

# plot
plot(x=res_shrink$log2FoldChange,
     y=log_padj,
     pch=20,
     cex=.5,
     col=(log_padj > 10)+1, # color padj < 0.1 red
     ylab="negative log-scaled adjusted p-value",
     xlab="shrunken log2 fold changes")

#############

# PCA plot

# normalized, variance-stabilized transformed counts for visualization
vsd <- vst(dds, blind=FALSE)

PC_Data<- plotPCA(vsd, intgroup=c("age"), returnData=TRUE)
PC_Data

rld <- rlogTransformation(dds, blind=T)
rv <- rowVars(assay(rld))
select <- order(rv, decreasing=T)[seq_len(min(500,length(rv)))]
pc <- prcomp(t(assay(vsd)[select,]))
pc
# alternatively, using ggplot

data <- plotPCA(vsd,returnData=TRUE,intgroup=c("age"))
data
p <- ggplot(data,aes(x=PC1,y=PC2,col=age))
p <- p + geom_point() + 
  xlab(paste("PC1: ", round(attr(data,"percentVar")[1],2)*100, "% variation explained", sep="")) + 
  ylab(paste("PC2: ", round(attr(data,"percentVar")[2],2)*100, "% variation explained", sep="")) +
  geom_label_repel(aes(label=name))
p


#write csv to get corresponding gene name from ensembl
# drop <- c("baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")
# ensemblnames <- results[,!(names(results) %in% drop)]
# ensemblnames <- results[0]
# ensemblnames
# write.csv(ensemblnames, "ensemblnames.csv")



#heatmap

# regularized log transformation of counts
  #rld <- rlog(deseq, blind=FALSE)

# order gene names by absolute value of shrunken log2 fold change (excluding cook's cutoff outliers)
lfcorder <- data.frame(res_shrink) %>%
  filter(!is.na(padj) & padj<0.2) %>% 
  arrange(-abs(log2FoldChange)) %>% 
  rownames() 

#ensemblnames <- lfcorder[1:30]
#write.table(ensemblnames, file="ensemblnames.txt", sep = "\t",
#            row.names = FALSE, col.names = FALSE, quote=FALSE)

humangenes <- read.csv("Human_reference_May1.csv")

# create a metadata data frame to add to the heatmaps
df <- data.frame(colData(dds)[,c("age")])
rownames(df) <- colnames(dds)
colnames(df) <- c("age")

#df_humangenes <- merge(x = df, y = humangenes, by.x= )
# use regularized log-scaled counts
#heatmapdf <- assay(rld)
heatmapdf <- as.data.frame(counts(dds,normalized =TRUE))

merged_data <- merge(x=heatmapdf,y=humangenes, by.x=0, by.y="GeneID", all.x=TRUE)

heatmap_metrics <- merged_data[merged_data$Row.names %in% lfcorder[1:30],]

#heatmap_metrics_merged <- merge(x=heatmap_metrics,y=humangenes, by.x=0, by.y="GeneID", all.x=TRUE)


pheatmap(
  as.matrix(heatmap_metrics[,2:10]), 
  cluster_rows=TRUE, 
  show_rownames=TRUE,
  cluster_cols=TRUE,
  labels_row=heatmap_metrics$Gene_Name,
  annotation_col=df,
  scale="row"
  )

# re-scale regularized log-scaled counts by baseMean (estimated mean across all samples)
pheatmap(
  assay(rld)[lfcorder[1:30],] - log(results[lfcorder[1:30],"baseMean"],2), 
  cluster_rows=TRUE, 
  show_rownames=TRUE,
  cluster_cols=TRUE,
  labels_row=
  annotation_col=df
)

# re-scale regularized log-scaled counts by reference level
pheatmap(
  assay(rld)[lfcorder[1:30],] - rowMeans(assay(rld)[lfcorder[1:30],dds$age=="young"]), 
  cluster_rows=TRUE, 
  show_rownames=TRUE,
  cluster_cols=TRUE,
  annotation_col=df
)

##vijender chunk of code

normcounts<-as.data.frame(counts(dds,normalized =TRUE))

normAnnot<-merge(normcounts,ref, by=0, all.x=T)
normAnnot = data.frame(lapply(normAnnot, as.character), stringsAsFactors=FALSE)
write.csv(normAnnot,file="NormalisedCounts.csv")