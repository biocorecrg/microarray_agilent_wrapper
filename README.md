# Microarray Agilent wrapper

## Goal:
Expression microarrays aim to study expression profiles of genes across different experimental groups, or over time.

## Tools:
Here, all analysis are performed in the R statistical environment and using Bioconductor packages (particularly the limma package).

## Wrapper:
The wrapper function will perform quality control and differential expression of a given Agilent microarray dataset.

### source file containing functions
```r
source("./wrapper_microarray_agilent.R")
```

function <b>agilent_wrap</b> takes the following arguments:
* Name given to the project
* Path to folder containing the raw data
* Path to the folder where you want the analysis to be saved
* target file, in a R object
* comparisons to be performed, in a R vector

and is called the following way:

```r
agilent_wrap("MyAnalysis", "/tmp/raw", "/tmp/analysis", target_file, comparisons)
```

### preparing the target file
The target file relates raw data file names to experimental groups.
It is made of 3 columns (2 mandatory and 1 optional): name is important and case sensitive, but order doesn't matter:
+ <b>"FileName"</b> contains the name (not the path!!) of the raw data files: one line represents one sample and one file.
+ <b>"Group"</b> contains the name of the experimental groups to which each sample belongs.
+ <b>"Sample"</b> (optional) contains the simplified and more comprehensive name of samples. Defaults to FileName.
+ <b>"Pairs"</b> (optional) contains the sample pairs (same animal, same patient, same plate etc.). If that column is present, a paired moderated t-test will be applied to account for that bias. No default.
+ <b>"Batch"</b> (optional) contains the batch number: can be the slide number or any technical parameter that you think can influence the day (day of the experiment etc.). If that column is present, an additional clustering plot will be produced that colors sample/file names per batch. No default.

### Example run
<b>prepare data</b>

```r
  # folder containing raw data (.txt files from Agilent/Feature extraction)
rawd <- "/tmp/data/"
	# folder where analysis will be written
outd <- "/tmp/analysis/"
	# sample names to input in target file
sampn <- substr(dir(rawd), 33, c(rep(40, 4), rep(39,4)))
	# group names to input in target file
groupn <- substr(dir(rawd), 33, 38)
	# target file
targ <- data.frame(FileName=dir(rawd), Group=groupn, Sample=sampn)
	# compute all possible comparisons: the command below produces a vector in the correct format
compall <- apply(combn(unique(groupn), 2), 2, function(x)paste(x, collapse="-"))
```

Note that if you prepare the target file in a text editor and not using regular expression, you can read it in the following way:
```r
targ <- read.table("targets.txt", sep="\t", header=T, as.is=T)
```

If you prepare the target file in Excel, the easiest is to save it in csv (comma separated) format and read it in the following way:
```r
targ <- read.csv("targets.csv", sep=",", header=T, as.is=T)
```

<b>launch analysis</b>
```r
agilent_wrap("TestAnalysis", rawd, outd, targ, compall)
```

The results are organised the following way:<br>

TestAnalysis/<br>
├── DE<br>
│   ├── DE_analysis.xls # pairwise differential expression analysis table: fold changes, p-values, annotation etc. <br>
│   └── DE_summary.txt # summary of differential expression analysis given changing criteria of fold change and adjusted p-values<br>
├── QC<br>
│   ├── Boxplots.pdf # boxplots of raw and normalized data<br>
│   ├── Clustering.pdf # hierarchical clustering of samples colored by sample group<br>
│   └── MAplots # MA plots<br>
└── target_file.txt # target file<br>

## Steps:
* Reading in data
* Data preprocessing:
** Background correction
** Normalization
* Quality Control of raw and normalized data
* Differential expression analysis

## Analysis of Agilent single channel gene expression arrays

### Reading in data
Raw data of Agilent microarrays consists of text files output of Feature Extraction software (Agilent software).
limma Bioconductor package provides a function for reading in intensities and annotation from these files.
<br>
<b>Our example set</b>
6 samples: 3 KO and 3 WT.

```r
# Create a target file: relate files to experimental groups
## column "FileName" contains the name fo the raw data files: one line represents one sample and one file.
## column "Cy3" contains the name of the experimental groups to which each sample belongs.

# target file can be simply build like:
targets <- data.frame(FileName=c("KO1_rawdata.txt", "KO2_rawdata.txt", "KO3_rawdata.txt", "WT1_rawdata.txt", "WT2_rawdata.txt", "WT3_rawdata.txt"), Cy3=c(rep("KO", 3),rep("WT", 3)))

# if you want to retrieve files names from a specific directory, and use regular expressions to retrieve experimental group names:
raw.files <- dir(path="RawDataDir", pattern=".txt$")
exp.groups <- gsub("[0-9]{1}_rawdata.txt", "", raw.files)

targets <- data.frame(FileName=raw.files, Cy3=exp.groups)

# Load limma library
library("limma")

# Read in the raw data in an object (here Eraw):
Eraw <- read.maimages(files=targets$FileName, source="agilent", path="RawDataDir", green.only=T, 
                      column=c(E="gMeanSignal", Eb="gBGMedianSignal"), 
                      annotation=c("FeatureNum", "Sequence", "ControlType", "ProbeName", "GeneName", "SystematicName", "Description"))

# E is the column that is read as the foreground signal.
# Eb is the column that is read as the background signal.

Read RawDataDir/KO1_rawdata.txt 
Read RawDataDir/KO2_rawdata.txt 
Read RawDataDir/KO3_rawdata.txt 
Read RawDataDir/WT1_rawdata.txt 
Read RawDataDir/WT2_rawdata.txt 
Read RawDataDir/WT3_rawdata.txt 

# The object obtained is of class "EListRaw".

Eraw$E # Foreground signal.

Eraw$Eb # Background signal.

Eraw$genes # Annotations
```

### Data preprocessing

<b>Background correction: normexp</b>

```r
E.bg <- backgroundCorrect(Eraw, method="normexp", offset=50)
Array 1 corrected
Array 2 corrected
Array 3 corrected
Array 4 corrected
Array 5 corrected
Array 6 corrected

```

<b>Quantile normalization</b>

```r
# performed on log2 transformed data
E.norm <- normalizeBetweenArrays(log2(E.bg$E), method="quantile")

# Modify column names of the E.norm object, for further usage (for example, plots).
colnames(E.norm) <- gsub("_rawdata", "", colnames(E.norm))

# Change rownames for the FeatureNum IDs (unique)
rownames(E.norm) <- Eraw$genes$FeatureNum
```

### Quality control of raw and normalized data

<b>Boxplots on raw and normalized data</b>

```r
pdf("Boxplots.pdf")
statboxRaw <- boxplot.stats(log2(Eraw$E), coef=3,do.conf = F, do.out = F)
statboxNorm <- boxplot.stats(E.norm, coef=3, do.conf = F, do.out = F)
boxplot(data.frame(log2(Eraw$E)), las = 2, main="Raw data", ylim=c(statboxRaw$stats[1],statboxRaw$stats[5]), outline=F, names=targets$Cy3)
boxplot(data.frame(E.norm), las = 2, main="Quantile Normalized data", ylim=c(statboxNorm$stats[1],statboxNorm$stats[5]), outline=F, names=targets$Cy3)
dev.off()
```

<br>

<b>MA-plots</b>
In wrapper
<br>

<b>Hierarchical clustering on normalized data</b>

```r
# Load ggdendro library
library(ggdendro)

# Perform the clustering on the most dynamic genes (coefficient of variation > 0.1)
std <- apply(E.norm, 1, sd) 
m <- rowMeans(E.norm)
use <- std / m > 0.1

# Build hierarchical clustering object (describes the tree produced by the clustering process, given the method applied)
hc <- hclust(d=dist(t(E.norm[use,])), method="ave") 

# save in pdf file (alternatively jpeg, png, tiff,...)
pdf("Clustering.pdf")
ggdendrogram(hc, rotate = FALSE, size = 3, leaf_labels=FALSE) + labs(title="Hierarchical clustering of normalized data")
dev.off()
```

<br>

### Differential expression analysis

```
# create design matrix
design <- model.matrix(~0+ targets$Cy3)
colnames(design) <- substr(colnames(design), 12, nchar(colnames(design)))

# Remove control probes for the analysis:
## 0: non-control probes
## 1: positive controls
## -1: negative controls
E.nc <- E.norm[Eraw$genes$ControlType==0,]

# linear model fit
fit <- lmFit(E.nc, design)

# contrast / comparison to perform
contrast.matrix <- makeContrasts(KO-WT, levels=design)

# fit contrasts
fit2 <- contrasts.fit(fit, contrast.matrix)

# compute moderated t-statistics, moderated F-statistic, and log-odds of differential expression by empirical Bayes moderation of the standard errors towards a common value
fit2 <- eBayes(fit2)

# extract a table of differential expression genes from a linear model fit.
top <- topTable(fit2, coef=1, adjust="fdr", number=nrow(fit$coefficients), sort.by="none")

        logFC  AveExpr          t    P.Value adj.P.Val         B
4 -0.76336788 8.826017 -2.9667070 0.02594329 0.1238667 -4.018685
5 -0.11290779 6.268845 -0.9375064 0.38569134 0.6248556 -6.652367
6  0.01171608 6.134226  0.1354357 0.89682215 0.9534825 -7.118002
7  0.03797648 6.248262  0.4045178 0.70025339 0.8482111 -7.034735
8  0.05668749 6.345256  0.4577584 0.66368917 0.8264726 -7.008847

# calculate actual fold change from log2ratio information (logFC)
lr <- top$logFC
fc <- sign(lr)*2^abs(lr)

# add annotation to the table
top.1 <- data.frame(FeatureNum=rownames(top), top)
top.all <- merge(top.1, as.data.frame(Eraw$genes), by="FeatureNum", all=F)

# write in Excel file
library("WriteXLS")
WriteXLS("top.all", ExcelFileName="./DE_analysis.xls", col.names=T, row.names=F, AdjWidth=F, BoldHeaderRow=T, FreezeRow=1, FreezeCol=1)
```


<b> Multiple comparisons</b>

```r
# create design matrix
design <- model.matrix(~0+ targets$Cy3)
colnames(design) <- substr(colnames(design), 12, nchar(colnames(design)))
 
# Prepare vector of contrasts to be calculated
contrast.names <- c("KO-WT", "WT-KO")

## Function

# run as: diffexp(Eraw, E.norm, design, contrast.names)

diffexp <- function(E_raw, E_normalized, design, contrasts_vector){

	# preparing data frame
	top.all <- E_raw$genes[E_raw$genes$ControlType==0,]

	# keep only non control probes
	En <- E_normalized[E_raw$genes$ControlType==0,]

	# linear model fit (as before)
	fit <- lmFit(En, design)

	for(i in 1:length(contrasts_vector)){
      ci <- contrasts_vector[i]
      contrast.matrix <- makeContrasts(contrasts=ci, levels=design)
      fit2 <- contrasts.fit(fit, contrast.matrix)
      fit2 <- eBayes(fit2)
      top <- topTable(fit2, coef=1, adjust="fdr", number=nrow(fit$coefficients), sort.by="none")
      val <- grep("Val", colnames(top))
      colnames(top)[val] <- paste(colnames(top)[val], gsub("-","_vs_",ci), sep="_")
      val <- grep("Val", colnames(top), value=T)
      lr <- grep("log", colnames(top))
      colnames(top)[lr] <- paste("log2ratio", gsub("-","_vs_",ci), sep="_")
      lr <- grep("log", colnames(top), value=T)
      fc <- paste("FC", gsub("-","_vs_",ci), sep="_")
      top[,fc] <- sign(top[,lr])*2^abs(top[,lr])
      top.all <- data.frame(top.all, top[,c(lr,fc,val)])
	}
	# write in Excel file
	library("WriteXLS")
	WriteXLS("top.all", ExcelFileName="./DE_analysis.xls", col.names=T, row.names=F, AdjWidth=F, BoldHeaderRow=T, FreezeRow=1, FreezeCol=1)
}
```

<b>Note on paired data analysis</b>
<br>
By default, limma performs an unpaired moderated t-test.
<br>
But samples can be paired.
<br>
For example, in the case of tumour samples vs healthy samples, it is possible that sample Tumour1 was extracted from the same patient as sample Healthy1, in which case they are naturally paired.
<br>
If you want to take into account the paired status of your samples in the differential expression analysis, you will need to modify both target file and design.

```r
# add a column "Pairs" to the targets file:
targets <- data.frame(FileName=raw.files, Cy3=exp.groups, Pairs=rep(c(1,2,3),2))

# paired model matrix
design <- model.matrix(~0+targets$Pairs+targets$Cy3)
colnames(design) <- gsub("targets\\$|Cy3", "", colnames(design))
```

And proceed as for the unpaired differential expression analysis.

### Batch effect correction

Correction for batch effects can be attempted using the [http://www.bu.edu/jlab/wp-assets/ComBat/Abstract.html ComBat] algorithm, implemented in the "sva" Bioconductor package.

## Analysis of Affymetrix single channel gene expression arrays

### Reading in data
Raw data of Affymetrix microarrays consists of "CEL" binary files.
affy Bioconductor package provides functions to read and preprocess Affymetrix data.

<br>
<b>Our example set</b>
6 samples: 3 KO and 3 WT.

```r
# load affy package
library(affy)

# read data in working directory </b>
Data <- ReadAffy(celfile.path="./RawDataDir")

# You obtain an AffyBatch object
```

### Preprocessing

Preprocessing of Affymetrix raw data is usually done using the <b>RMA (Robust Multichip Average) method</b>, that consists of three steps:
* a background adjustment
* quantile normalization (see the Bolstad et al reference)
* summarization of probe sets

```r
# RMA preprocessing
eset <- rma(Data)

# extract raw data
raw.data <- exprs(Data)

# extract normalized data
affy.norm <- exprs(eset)
```

#### Quality control of raw and normalized data

<b>Histogram of intensities </b>

```r
hist(log2(raw.data))
hist(affy.norm)
```

<b>Affymetrix RNA degradation plot </b>
<br>
Used to detect possible RNA degradation.

```r
deg <- AffyRNAdeg(Data)
summaryAffyRNAdeg(deg)

# plot degradation
pdf("RNAdegradation.pdf")
plotAffyRNAdeg(deg)
dev.off()
```

<b> MA plots </b>
<br>
affy package contains a built in function to produce MA plots.

```r
jpeg("MAplots.jpg", width = 1960, height = 1960)
MAplot(Data,cex=3, pairs = TRUE)
dev.off()
```
<b> Boxplot of intensities</b>
<br>
<i>Same as for Agilent data</i>

<b> Hierarchical clustering </b>
<br>
<i>Same as for Agilent data</i>

#### Differential expression analysis
<i> Same as for Agilent data.</i>
