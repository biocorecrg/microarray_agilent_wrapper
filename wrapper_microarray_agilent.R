#### Wrapper function for microarray analysis of Agilent 1 color arrays
	# run as agilent_wrap("test", rawdir, outdir, targets, c("KO-WT"))

agilent_wrap <- function(namexp, rawpath, outpath, targetfile, comps){
	# check if needed packages are installed: if not, install them.
	list.of.packages <- c("ggplot2", "ggdendro","limma", "WriteXLS", "affy")
	new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
	if(length(new.packages)>0){install.packages(new.packages)}
	# load required packages
	require("limma")
	require("ggdendro")
	require("ggplot2")
	require("cowplot")
	require("affy")
	# get current directory (we will get back to it when the analysis is over)
	currpath <- getwd()
	# set output directory
	setwd(outpath)
	# create subdirectories
	dir.create(file.path(namexp), showWarnings = FALSE)
	dir.create(file.path(namexp, "QC"), showWarnings = FALSE)
	dir.create(file.path(namexp, "DE"), showWarnings = FALSE)
	dir.create(file.path(paste(namexp, "QC", sep="/"), "MAplots"), showWarnings = FALSE)
	# Read in the raw data in an object (here Eraw):
	Eraw <- read.maimages(files=targetfile$FileName, source="agilent", path=rawpath, green.only=T, column=c(E="gMeanSignal", Eb="gBGMedianSignal"), annotation=c("FeatureNum", "Sequence", "ControlType", "ProbeName", "GeneName", "SystematicName", "Description"))
	# Background correction: normexp 
	E.bg <- backgroundCorrect(Eraw, method="normexp", offset=50)
	# Quantile normalization performed on log2 transformed data
	E.norm <- normalizeBetweenArrays(log2(E.bg$E), method="quantile")
	rownames(E.norm) <- Eraw$genes$FeatureNum
	# If "targets" contains a column with sample name (additionally to file name), change raw and normalized data column names (also if sample names are unique)
	if(!is.null(targetfile$Sample) & length(unique(targetfile$Sample))==length(unique(colnames(E.norm)))){
		colnames(E.norm) <- targetfile$Sample
		colnames(Eraw$E) <- targetfile$Sample	
	}
	# save normalized intensities in file
	E.norm2 <- data.frame(Eraw$genes, E.norm)
	write.table(E.norm2, paste0(namexp, "/",Sys.Date(),"_",namexp,"_Normalized_intensities.txt"),sep="\t", col.names=T, row.names=F, quote=F)
	require("WriteXLS")
	WriteXLS("E.norm2", ExcelFileName=paste0(namexp, "/",Sys.Date(),"_",namexp,"_Normalized_intensities.xls"), col.names=T, row.names=F, AdjWidth=F, BoldHeaderRow=T, FreezeRow=1, FreezeCol=1)
	# set width of pdf plots given the number of samples...
	nsamp <- ncol(E.norm)
	pdfwd <- ifelse(nsamp <= 7, 7, 
			ifelse(nsamp > 7 & nsamp <= 12, 10, 
			ifelse(nsamp > 12 & nsamp <= 20, 14, 
			ifelse(nsamp > 20 & nsamp <= 30, 20, 
			ifelse(nsamp > 30 & nsamp <= 40, 26, 
			ifelse(nsamp > 40 & nsamp <= 50, 35, 40))))))
	# Boxplot of raw and normalized data
	pdf(paste0(namexp, "/QC/",Sys.Date(),"_",namexp,"_Boxplots.pdf"), width=pdfwd)
	statboxRaw <- boxplot.stats(log2(Eraw$E), coef=3,do.conf = F, do.out = F)
	statboxNorm <- boxplot.stats(E.norm, coef=3, do.conf = F, do.out = F)
	par(oma=c(5,2,2,2))
	boxplot(data.frame(log2(Eraw$E)), las = 2, main="Raw data", ylim=c(statboxRaw$stats[1],statboxRaw$stats[5]), outline=F)
	boxplot(data.frame(E.norm), las = 2, main="Quantile Normalized data", ylim=c(statboxNorm$stats	[1],statboxNorm$stats[5]), outline=F)
	dev.off()
	# Clustering
	# Perform the clustering on the most dynamic genes (coefficient of variation > 0.1)
	std <- apply(E.norm, 1, sd) 
	m <- rowMeans(E.norm)
	use <- std / m > 0.05
	# Build hierarchical clustering object (describes the tree produced by the clustering process, 	given the method applied)
	#hc <- hclust(d=dist(t(E.norm[use,])), method="ave") 
	# save in pdf file
	#cl <- ggdendrogram(hc, rotate = FALSE, size = 3, leaf_labels=FALSE) + labs(title="Hierarchical 	clustering of normalized data")

	dd.row <- as.dendrogram(hclust(dist(t(E.norm[use,]))), method="ave")
	ddata_x <- dendro_data(dd.row)
	hclab <- label(ddata_x)
	hclab$group <- NA
	for(g in unique(targetfile$Group)){
		if(!is.null(targetfile$Sample) & length(unique(targetfile$Sample))==length(unique(colnames(E.norm)))){
			gtmp <- targetfile$Sample[targetfile$Group==g]
		}else{
			gtmp <- targetfile$FileName[targetfile$Sample==g]
		}
		hclab$group[hclab$label%in%gtmp] <- g
	}
	cl <- ggplot(segment(ddata_x)) + geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) + geom_text(data=label(ddata_x), aes(label=label, x=x, y=0, colour=hclab$group), angle=90, hjust=1, size=4) + scale_y_continuous(limits=c(-max(segment(ddata_x)$y)/4, (max(segment(ddata_x)$y)+0.05*(max(segment(ddata_x)$y))))) + ggtitle(paste("Clustering of", namexp, "microarray project")) + scale_colour_discrete(name = "Groups")
	pdf(paste0(namexp, "/QC/",Sys.Date(),"_",namexp,"_Clustering.pdf"), width=pdfwd)
		plot(cl)
	dev.off()
	# dendrogram with colored batches, if "Batch" column is present
	if("Batch" %in% colnames(targetfile)){
		if( "Sample" %in% colnames(targetfile)){	
		hclab2 <- merge(hclab, targetfile[,grep("Sample|Batch", colnames(targetfile))], by.x="label", by.y="Sample")
 		hclab <- hclab2[order(hclab2$x),c("x", "y", "label", "group", "Batch")]
		}else{
		hclab2 <- merge(hclab, targetfile[,grep("FileName|Batch", colnames(targetfile))], by.x="label", by.y="FileName")
 		hclab <- hclab2[order(hclab2$x),c("x", "y", "label", "group", "Batch")]
		}
	cl <- ggplot(segment(ddata_x)) + geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) + geom_text(data=label(ddata_x), aes(label=label, x=x, y=0, colour=hclab$Batch), angle=90, hjust=1, size=4) + scale_y_continuous(limits=c(-max(segment(ddata_x)$y)/4, (max(segment(ddata_x)$y)+0.05*(max(segment(ddata_x)$y))))) + ggtitle(paste("Clustering of", namexp, "microarray project colored by batch")) + scale_colour_discrete(name = "Batches")
	pdf(paste(namexp, "QC/Clustering_colbatches.pdf", sep="/"), width=pdfwd)
		plot(cl)
	dev.off()

	}

	# MA plots
	maplotfun(Eraw, E.norm, namexp)
	# differential expression analysis: call diffexp function


	if("Pairs" %in% colnames(targetfile))
	{
  	cat("Paired analysis\n");
	# paired model matrix
	design <- model.matrix(~0+targetfile$Group+targetfile$Pairs)
	}else{
	design <- model.matrix(~0+ targetfile$Group)
	}
	colnames(design) <- gsub("targetfile\\$Group", "", colnames(design))
	colnames(design) <- gsub("targetfile\\$Pairs", "p", colnames(design))
	diffexp(Eraw, E.norm, design, comps, namexp)
	# save target file for potential later usage
	write.table(targetfile, paste0(namexp, "/",Sys.Date(),"_",namexp,"target_file.txt"), sep="\t", quote=F, row.name=F, col.name=T)
	# go back to previous directory
	setwd(currpath)
	
}

 
## Differential expression analysis function
	# run as: diffexp(Eraw, E.norm, design, contrast.names)
 
diffexp <- function(E_raw, E_normalized, design, contrasts_vector, name_exp){
	# preparing data frame
	top.all <- E_raw$genes[E_raw$genes$ControlType==0,]
	# keep only non control probes
	En <- E_normalized[E_raw$genes$ControlType==0,]
	# linear model fit (as before)
	fit <- lmFit(En, design)
 	# loop around contrasts
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
	# DE summary
	diffexpsum(top.all, name_exp, contrasts_vector)
	# write in Excel file
	require("WriteXLS")
	WriteXLS("top.all", ExcelFileName=paste0(name_exp, "/DE/",Sys.Date(),"_",name_exp,"DE_analysis.xls"), col.names=T, row.names=F, AdjWidth=F, BoldHeaderRow=T, FreezeRow=1, FreezeCol=1)
	# also write in text file for further manipulation of the data
	write.table(top.all, paste0(name_exp, "/DE/",Sys.Date(),"_",name_exp,"DE_analysis.txt"), quote=F, sep="\t", row.names=F, col.names=T)
}


## Differential expression summary function from limma's toptable
	# run as diffexpsum(toplimma, namexp, comparisons)
diffexpsum <- function(dataDif, name_exp, contrast){  
  fc <- grep("^FC", colnames(dataDif), value=T)
  qval <- grep("^adj", colnames(dataDif), value=T)
  
  contrast.names <- unlist(strsplit(contrast," {0,1}[,;] {0,1}"))
   
  diff.exp <- matrix(NA, ncol=length(contrast.names), nrow=8)

  lpval <- c(0.05, 0.01)
  lfc <- c(1.2, 2, 5, 10)
  cnames <- list()
 
  for(cont in 1:length(contrast.names)){
  k <- 1
      for(j in 1:length(lfc)){
        for(i in 1:length(lpval)){
        diff.exp[k,cont] <- length(which(abs(dataDif[,fc[cont]]) > lfc[j] & dataDif[,qval[cont]] < lpval[i]))
        cnames[k] <- paste("|FC| > ", lfc[j], " AND adjPval < ", lpval[i], sep="") 
        k <- k+1
      }
    }
  }
  
  colnames(diff.exp) <- contrast.names
  rownames(diff.exp) <- cnames
  write.table(diff.exp, paste0(name_exp, "/DE/",Sys.Date(),"_",name_exp,"DE_summary.txt"), sep="\t", col.names=NA)
}


## MA plots function
	# run as maplotfunc(RawDataObject, NormDataObject)
maplotfun <- function(DataRaw, DataNorm, name_exp){
	MAfun <- function(data, j, gplot, log=FALSE, ...){
  		if(log){
    			data <- log2(data)
  		}
  		# 'median array'
  		med <- apply(data, 1, median)
  		# M and A values
  		M <- data[,j] - med
  		A <- (data[,j] + med) / 2
  		gmain <- paste(colnames(data)[j],gplot,sep=" ")
  		ma.plot(A, M, main=gmain, cex=1, pch=19)
		}
	# before normalization
	tmp <- DataRaw$E
	colnames(tmp) <- gsub(".txt","",colnames(tmp))
	lf <- length(colnames(tmp))
	for(i in 1:lf){
  		ma_af <- paste(name_exp, "/QC/MAplots/MAplots_beforeNormalization_", colnames(tmp)[i], ".png",sep="")
  		png(ma_af, width = 360, height = 360)
  		MAfun(tmp, i, "Before Normalization",log=T, pch=16, show.statistics=F)
   		dev.off()
	}                                                                                           
	# after normalization
	tmp <- DataNorm
	colnames(tmp) <- gsub(".txt","",colnames(tmp))
	for(i in 1:lf){
  		ma_bef <- paste(name_exp, "/QC/MAplots/MAplots_afterNormalization_", colnames(tmp)[i], ".png",sep="")
  		png(ma_bef, width = 360, height = 360)
  		MAfun(tmp, i, "After Normalization", pch=16, show.statistics=F)
  	dev.off()
	}
}





