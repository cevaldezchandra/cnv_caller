##### cfDNA CN caller - ########################
### Written by Crystal Valdez 
### Date: August 15, 2017
### Version: 2.0
###############################################

### Setting input paths for EXCAVATOR calculations
# Input should look like:
# Rscript CNV_v7.R /path/to/sample.bam targetname /path/to/output
#
# Takes in 3 commandline prompts: 
# 1. Path to .bam file
# 2. Target name used to name target initialization folder and output folder
# 3. Path to desired output directory - WITHOUT folder name, the script will
#    generate individual folders using the targetname

### Input files and data needed (bullets 2-4 are found in main EXCAVATOR folder)
# 1. Sample.RCNorm_txt (output from EXCAVATOR)
# 2. 51_total_exons.csv (51 unpaired samples for pooled background)
# 3. gene_names.txt
# 4. exon_names.txt

### Output files produced


vars <- commandArgs(trailingOnly = TRUE)

### Setting input paths for EXCAVATOR calculations
BamFile <- vars[1]
TargetName <- vars[2]
Output <- vars[3] # with batch#
#Cnv_output <- vars[4] # added this for cnv analysis
Background <- vars[4] # choose background to use
RunExac <- vars[5] # flag is EXCAVATOR should be run

### Create necessary files, folders and paths
PathEXAC <- file.path("/usr/applinux/tools/EXCAVATOR_Package_v2.2/EXCAVATOR")
dir.create(file.path(Output,TargetName), showWarnings = FALSE)
setwd(file.path(Output, TargetName)) # set output as wd - NEED both lines
wd <- setwd(file.path(Output, TargetName)) #set output as working directory

# copy over exacavtor directory
#system(paste("cp -r", PathEXAC, wd))
#print(paste("Working directory: ",wd))

# copy over exacavtor directory for calculations, if needed
#ifelse(!dir.exists(file.path(mainDir, subDir)), dir.create(file.path(mainDir, subDir)), FALSE)
ifelse(dir.create(file.path(wd, "EXCAVATOR")), system(paste("cp -r", PathEXAC, wd)), FALSE)
#ifelse(!dir.exists(file.path(wd, "EXCAVATOR")), system(paste("cp -r", PathEXAC, wd)), FALSE)
PathEXAC2 <- file.path(wd, "EXCAVATOR")
print(paste("New EXCAVATOR path: ", PathEXAC2))

### Permanant place for EXACAVTOR program folders, .bed file and SourceTarget.txt
BedFile <- file.path(PathEXAC2,"cfDNA_panel.bed")
SourceFile <- file.path(PathEXAC2,"SourceTarget.txt")
Target <- file.path(PathEXAC2,"TargetPerla.pl")
Read <- file.path(PathEXAC2,"ReadPerla.pl")
ReadInput <- cat(c(TargetName, "hg19", BamFile, TargetName, "CT"), file="ReadInput.txt", sep=" ")
ReadInput <- file.path(wd,"ReadInput.txt")
PooledBack <- file.path(PathEXAC2, "51_total_exons.csv")
print(paste("BELOW ReadInput.txt",ReadInput))
#Output <- file.path(pwd) # for running locally on machine
print(paste(".bam file path: ", BamFile))
print(paste("bed file path: ", BedFile))
print(paste("SourceTarget.txt file located: ", SourceFile))
print(paste("Sample name: ", TargetName))
print(paste("Output folder: ", Output,"/",TargetName,sep=""))

# new output folders
#dir.create(file.path(Cnv_output, TargetName)) # new cnv output
#setwd(file.path(Cnv_output, TargetName))
#out <- setwd(file.path(Cnv_output, TargetName))

#ref_pooled <- "/work/Paloma_CF/cnv_51_unpaired/ref_pooled" # hardcoded to include various pooled means

print(paste("Sample name: ", TargetName))
print(paste("EXCAVATOR folder: ", Output,"/",TargetName,sep=""))
print(paste("CNV pooled output folder: ", wd))
print(paste("Current working directory: ",getwd()))

### Libraries needed by R: need Hmisc R package for EXCAVATOR
list.of.packages <- c("Hmisc")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos='https://cran.r-project.org')
library(Hmisc)

### Run EXCAVATOR  perl scripts
cmd1 <- paste(Target, SourceFile, BedFile, TargetName)
cmd2 <- paste(Read, ReadInput, wd, "--mode pooling")

print(paste("EXCAVATOR step 1: ",cmd1))
print(paste("EXCAVATOR step 2: ",cmd2))
# check is EXCAVATOR has been run already or not
if (RunExac == "T") {
    system(cmd1)
    system(cmd2)
    print("Running EXCAVATOR - part I")
}  else {
    print("EXCAVATOR has already run, continue...")
}
print(paste("directory where normalization of RC took place: ", getwd()))
print("EXCAVATOR RCNorm generating step complete....")
###### EXCAVATOR COMPLETE#########


###### CNV caller ################
### Read in Output files from EXCAVATOR
samp <- list.files(path=wd, pattern=".RCNorm_txt", full.names = TRUE)
samples <- samp[!grepl("Control.RCNorm_txt",samp)] # print all but control
print(samples)

### Extract out only RCNorm column, field 6, and normalize samples across all exons 
RC <- lapply(samples, read.table, header=T, colClasses=c(rep("NULL",5),"numeric"), sep="\t")
Ex <- lapply(samples, read.table, header=T, colClasses=c(rep("NULL",4),"character","NULL"),sep="\t")
exon <- read.table(paste(PathEXAC2,"/exon_names.txt",sep=""), sep="\t", col.names = "exon")
xcn <- data.frame(RC) # need to convert to data.frame
df = as.matrix(xcn)
preN <- t(t(df)/colSums(df)) # divide by sum of column, do not multiply by 1mil here
norm <- data.frame(preN)
RCn <- norm*1000000 # multiply by 1mil, assume 1mil reads

### make exon names to data.frame rownames
NcnE <- cbind(Ex,RCn)
rownames(NcnE) <- NcnE[,1]
NcnE[,1] <- NULL
cn=NcnE

xcnE <- cbind(Ex,xcn)
rownames(xcnE) <- xcnE[,1]
xcnE[,1] <- NULL

write.table(xcn, file="RC.xls",sep="\t", row.names=TRUE, col.names=TRUE)
write.table(RCn, file="Normalized_RC.xls",sep="\t", row.names=TRUE, col.names=TRUE)

### Bring in background samples for calculating ratios

#back_14 <- read.table(paste(ref_pooled,"/",Background,sep=""),sep=",", row.names = 1, header = TRUE)
back_14 <- read.table(paste(PathEXAC2,"/",Background,sep=""),sep=",", row.names = 1, header = TRUE)
back_14_adj <- t(t(back_14)/colSums(back_14))*1000000

### remove exons that are not captures in the sample
back_clean <- back_14_adj[(rownames(back_14_adj) %in% rownames(cn)),]
back_14_adj <- back_clean

### Trimming the background mean by removing top and bottom 10% of exons 
# added by Crystal Valdez July 20, 2017
# changed by Crystal Valdez August 11, 2017 - taking 10% of 51 samples
f <- function(x) {
  or <- sort(x)
  trim_t <- tail(or, -5) #adjust according to number of exons to remove
  trim_ht <- head(trim_t, -5) #adjust according to number of exons to remove
  mean_ht <- mean(trim_ht)
}
trim_mean <- as.data.frame(apply(back_14_adj,1,f))

# calculate new standard deviation
trim_sd <- as.data.frame(apply(back_14_adj,1,function (x){
  or <- sort(x)
  trim_t <- tail(or, -5) #adjust as needed
  trim_ht <- head(trim_t, -5) #adjust as needed
  st_d <- sd(trim_ht)
}))


###### Ready for input sample
#med=apply(log2(back_med),2,median)
med=mean(apply(log2(trim_mean),2,median)) ###mean of medians from background

log.rT=cn
for(i in 1:ncol(cn)) {
  #log.rT[,i]=log2(cn[,i])-apply(log2(trim_mean),1,mean)+med-median(log2(cn[,i]))
  log.rT[,i]=log2(cn[,i])-apply(log2(trim_mean),1,mean) #removed median centering Crystal Valdez 17AUG2017
}

### Calculate Z-scores
log.rT_linear = 2^(log.rT)
zz=(log.rT_linear-trim_mean)/trim_sd

# added print all exons and Z-scores
write.table(log.rT, file="exons_log2ratio_all_exons.xls", sep="\t", row.names=TRUE)
write.table(zz, file="zscores_all_exons.xls", sep="\t", row.names=TRUE)
### calculate pnorm values for each exon old code
#loocv_norm2 <- apply(log.rT,1,function(x) pnorm(x,mean=0,sd=1))
#write.table(loocv_norm2, file="pnorm_all_exons.xls", sep="\t", row.names=TRUE) 

### Preparing for gene level summary 
# extract out exons with very low ANOVA scores - updated August 16, 2017
pvalue <- file.path(PathEXAC2, "p_value_remove.txt")
exonn <- read.table(pvalue, sep="\t", row.names = 1, header=TRUE)
#exon.55 <- data.frame(exonn[,1])
#exrow <- cbind(exon.55, exon.55) # need to create row.names
#rownames(exrow) <- exrow[,1]
#exrow[,1] <- NULL
#print(exrow)

### Remove roughly 10% (86) exons determined by ANOVA analysis to be too variable
rownames(exonn) -> remove
exon.c <- log.rT[!rownames(log.rT) %in% remove, , drop=FALSE] #remove 86 exons 
z.c <- zz[!rownames(zz) %in% remove, ,drop=FALSE] # remove 86 exons z-scores
#print(exon.c)
#boxplot(exon.c, las=2, main="Log2(Ratio Sample/Background) - 55 exons removed")
#stripchart(exon.c, vertical=TRUE, pch=21, col=c("skyblue","steelblue","skyblue4"),
#           bg="black", method="jitter", add=TRUE)
write.table(exon.c, file="exons_log2ratio_REMOVED.xls", sep="\t", row.names=TRUE) 


### Extracing gene info from cleaned up exons
genelist <- file(paste(PathEXAC2,"/gene_names.txt",sep=""),open="r")
linn <- readLines(genelist)

totalg = NULL # for mean - 55 removed
totals = NULL # for median - 55 removed
totalz = NULL # for Z-scores - 55 removed
Atotalg = NULL # for mean - all exons
Atotals = NULL # for median - all exons
Atotalz = NULL # for Z-scores - all exons
for (i in linn){
  # group exons according to gene name - 55 removed exons
  gene <- exon.c[grep(i, rownames(exon.c)), , drop=FALSE]
  zscor <- z.c[grep(i, rownames(z.c)), , drop=FALSE]
  # group exons accoding to gene name - all
  Agene <- log.rT[grep(i,rownames(log.rT)), , drop=FALSE]
  Azscor <- zz[grep(i, rownames(zz)), , drop=FALSE]
  # calculate mean, median
  geneT <- colMeans(gene)
  geneS <- apply(gene,2,median)
  zscorT <- colMeans(zscor)
  
  AgeneT <- colMeans(Agene)
  AgeneS <- apply(Agene,2,median)
  AzscorT <- colMeans(Azscor)
  # combind all genes together
  totalg <- rbind(totalg, geneT)
  totals <- rbind(totals, geneS)
  totalz <- rbind(totalz, zscorT)
  
  Atotalg <- rbind(Atotalg, AgeneT)
  Atotals <- rbind(Atotals, AgeneS)
  Atotalz <- rbind(Atotalz,AzscorT)
  do.call(rbind, exon.c)
}

### Boxplot of Gene Summaries and Z-scores
totalg <- data.frame(totalg)
boxplot(totalg, las=2, main="Gene Level Summary - REMOVED exons", outline=FALSE, ylim=c(-1,1))
stripchart(totalg, vertical=TRUE, pch=21, col=c("skyblue","steelblue","skyblue4"),
           bg="black", method="jitter", add=TRUE)
abline(h=0.2)
abline(h=-0.2)

Atotalg <- data.frame(Atotalg)
boxplot(Atotalg, las=2, main="Gene Level Summary - all exons", outline=FALSE, ylim=c(-1,1))
stripchart(Atotalg, vertical=TRUE, pch=21, col=c("skyblue","steelblue","skyblue4"),
           bg="black", method="jitter", add=TRUE)
abline(h=0.2)
abline(h=-0.2)

totalz <- data.frame(totalz)
boxplot(totalz, las=2, main="Z-score Summary/Gene - REMOVED exons", outline=FALSE, ylim=c(-1,1))
stripchart(totalz, vertical=TRUE, pch=21, col=c("skyblue","steelblue","skyblue4"),
           bg="black", method="jitter", add=TRUE)
abline(h=0.2)
abline(h=-0.2)

Atotalz <- data.frame(Atotalz)
boxplot(Atotalz, las=2, main="Z-score Summary/Gene - all exons", outline=FALSE, ylim=c(-1,1))
stripchart(Atotalz, vertical=TRUE, pch=21, col=c("skyblue","steelblue","skyblue4"),
           bg="black", method="jitter", add=TRUE)
abline(h=0.2)
abline(h=-0.2)

### Add gene names as rownames 
Gname <- cbind(linn,totalg)
totalG <- data.frame(Gname)
print(totalG)
rownames(totalG) <- totalG[,1]
totalG[,1] <- NULL
write.table(totalG, file="gene_summary_REMOVED.xls", sep="\t")

Gname <- cbind(linn,Atotalg)
AtotalG <- data.frame(Gname)
print(AtotalG)
rownames(AtotalG) <- AtotalG[,1]
AtotalG[,1] <- NULL
write.table(AtotalG, file="gene_summary_all_exons.xls", sep="\t")

### Add gene names to mean of Z-scores
Zname <- cbind(linn,totalz)
totalzz <- data.frame(Zname)
rownames(totalzz) <- totalzz[,1]
totalzz[,1] <- NULL
write.table(totalzz, file="zscore_gene_summary_REMOVED.xls", sep="\t")

Zname <- cbind(linn,Atotalz)
Atotalzz <- data.frame(Zname)
rownames(Atotalzz) <- Atotalzz[,1]
Atotalzz[,1] <- NULL
write.table(Atotalzz, file="zscore_gene_summary_all_exons.xls", sep="\t")


#write.table(totalGene, file="Gene_Summary_mean.xls", sep="\t", row.names=TRUE) # will need to add genes later
closeAllConnections()
print("Gene level summary complete")
#########################


