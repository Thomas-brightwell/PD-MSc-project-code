#Set wd and load in the required modules
setwd("C:/Users/Thomas Brightwell/OneDrive - Imperial College London/PD_MSc_Project2024/Working materials/Visual Studio/eQTL")

install.packages("tidyverse")
install.packages("MatrixEQTL")
install.packages("remotes")
remotes::install_github("matthiasheinig/eQTLpipeline")
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("VariantAnnotation")
BiocManager::install('snpStats')
install.packages("LDlinkR")
BiocManager::install("biomaRt")
install.packages("ActivePathways")
install.packages("genetics")
BiocManager::install("limma")
BiocManager::install("DESeq2")
library(snpStats)
library(tidyverse)
library(MatrixEQTL)
library(VariantAnnotation)
library(eQTLpipeline)
library(LDlinkR)
library(biomaRt)
library(ActivePathways)
library(genetics)
library(limma)
library(DESeq2)
#First,load in and cross-reference the lead SNPs, to see which are missing
leadSNPs <- read.table("TableS2_v2.txt", header = TRUE, sep = "\t")
lookuptable <- read.table("Lookuptable.txt", header = TRUE, sep = "\t")
colnames(leadSNPs)
colnames(lookuptable)
leadSNPsannotated <- inner_join(leadSNPs, lookuptable, join_by("snp" == "rs_id_dbSNP151_GRCh38p7"))
missingleadSNPs <- anti_join(leadSNPs, lookuptable, join_by("snp" == "rs_id_dbSNP151_GRCh38p7"))
test <- lookuptable[lookuptable$rs_id_dbSNP151_GRCh38p7 == "rs114138760",]
rm(leadSNPsannotated,missingleadSNPs, test)
#rs76763715 is matched with rs560496384
#rs5019538 is matched with rs356216
#rs112485576 is a synonym of rs504594, which is in the WGS data we have
#rs34637584 is a pathogenic variant, so we can't use it as the frequency is too low

#Now, add those SNPs into the table for use
leadSNPs <- subset(leadSNPs, select = "snp")
leadSNPs$snp <- sub("rs112485576", "rs504594", leadSNPs$snp)
leadSNPs$snp <- sub("rs76763715", "rs560496384", leadSNPs$snp)
leadSNPs$snp <- sub("rs5019538", "rs356216", leadSNPs$snp)
leadSNPs <- inner_join(leadSNPs, lookuptable, join_by("snp" == "rs_id_dbSNP151_GRCh38p7"))





#load the coordinates of the LD blocks and add the LD block number to the lead SNP table
blockcoordsb38 <- read.table("blockcoordsb38.bed", header = TRUE, sep = "\t")
colnames(blockcoordsb38) <- c("chr", "start", "end")
blockcoordsb38$LDblock <- 1:nrow(blockcoordsb38)
leadSNPs$LDblock <- c(1:6, 6:57, 59:79, 79:88)
leadSNPs2 <- inner_join(leadSNPs, blockcoordsb38, join_by("LDblock" == "LDblock", "chr" == "chr"))
leadSNPs2 <- leadSNPs2 %>%
    filter(variant_pos >= start & variant_pos <= end)

#Load in the genotype file for the loci
genotypeb38 <- read.table("PDlocib38.recode.vcf", sep = "\t", header = TRUE)
genotypeb38 <- inner_join(genotypeb38, lookuptable, join_by("ID" == "variant_id"))
genoinfo <- genotypeb38[, 1:9]
genoinfo <- cbind(genoinfo, genotypeb38[, 848:854])
#Classify if the variants are SNPs or not
genoinfo$ref_length <- nchar(genoinfo$ref)
genoinfo$alt_length <- nchar(genoinfo$alt)
genoinfo$SNP <- ifelse(genoinfo$ref_length == 1 & genoinfo$alt_length == 1, 1, 0)
table(genoinfo$SNP)
#Add a tag to the lead SNPs
genoinfo$lead <- ifelse(genoinfo$ID %in% leadSNPs$variant_id, 1, 0)
table(genoinfo$lead)
#Convert the genotype data to a 0 1 2 format
genoinfo <- cbind(genoinfo, genotypeb38[, 10:847])
for (column in names(genoinfo)[21:858]) {
    genoinfo[[column]] <- sub(":.*", "", genoinfo[[column]])
    genoinfo[[column]] <- ifelse(genoinfo[[column]] == "0/0", 0, ifelse(genoinfo[[column]] == "0/1", 1, 2))
}
#Add the LD block number to the genotype data
genoinfo$LDblock <- blockcoordsb38$LDblock[blockcoordsb38$chr == genoinfo$chr & genoinfo$variant_pos >= blockcoordsb38$start & genoinfo$variant_pos <= blockcoordsb38$end]
blocks <- left_join(genoinfo[,1:20], blockcoordsb38, by = "chr")
blocks <- blocks[blocks$variant_pos >= blocks$start & blocks$variant_pos <= blocks$end,]
genoinfo <- cbind(blocks, genoinfo[, 21:858])
table(genoinfo$LDblock)
table(blocks$lead)
write.table(genoinfo, "genotypeb38annotated.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

#Subset the genotype data to just europeans
genoinfo <- read.table("genotypeb38annotated.txt", sep = "\t", header = TRUE)
phenotypes <- read.csv("GTExsubjectphenotypes.csv", header = TRUE)
phenotypes <- phenotypes[phenotypes$RACE == 3,]
info <- genoinfo[, 1:23]
samples <- genoinfo[, which(colnames(genoinfo) %in% phenotypes$SUBJID)]
genoinfo <- cbind(info, samples)
rm(info, samples)

#calculate the alternate allele frequency
genoinfo$ALT_FREQ <- rowSums(genoinfo[, 24:ncol(genoinfo)])/(2*ncol(genoinfo[, 24:ncol(genoinfo)]))
ALTFREQ <- genoinfo$ALT_FREQ
genoinfo <- cbind(genoinfo[,1:23], ALTFREQ, genoinfo[,24:738])
rm(ALTFREQ)


summary_genoinfo <- genoinfo %>% group_by(LDblock) %>% summarise(count = n(), avefreq = mean(ALTFREQ), minfreq = min(ALTFREQ), maxfreq = max(ALTFREQ), quant0.01 = quantile(ALTFREQ, 0.01), quant0.99 = quantile(ALTFREQ, 0.99), zerocount = sum(ALTFREQ == 0))
write.csv(summary_genoinfo, "altfreqinfo.csv", row.names = FALSE)

#calculate the correlation between the lead SNP and the other SNPs in the LD block
genoinput <- genoinfo[, -which(names(genoinfo) %in% c("CHROM", "ALTFREQ", "POS", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "chr", "variant_pos", "ref", "alt", "num_alt_per_site", "rs_id_dbSNP151_GRCh38p7", "ref_length", "alt_length", "SNP", "start", "end", "variant_id_b37", "missing"))]
rownames(genoinput) <- genoinput$ID
genoinput <- genoinput[, 2:718]



#block 6 and 79 have 2 lead SNPs, so need to do these blocks separately. Block 58 is also separated due to the pathogenic lead SNP.
block6 <- genoinput[genoinput$LDblock == 6, 3:717]
block79 <- genoinput[genoinput$LDblock == 79, 3:717]
lead6p1 <- block6[which(rownames(block6) == "chr1_205754444_C_T_b38"),]
lead6p2 <- block6[which(rownames(block6) == "chr1_205768611_G_A_b38"),]
lead79p1 <- block79[which(rownames(block79) == "chr17_45666837_C_T_b38"),]
lead79p2 <- block79[which(rownames(block79) == "chr17_45720942_G_A_b38"),]
lead6p1 <- t(lead6p1)
lead6p2 <- t(lead6p2)
lead79p1 <- t(lead79p1)
lead79p2 <- t(lead79p2)
block6 <- t(block6)
block79 <- t(block79)
blockresults6p1 <- data.frame(2)
blockresults6p2 <- data.frame(2)
blockresults79p1 <- data.frame(2)
blockresults79p2 <- data.frame(2)
for (column in 1:ncol(block6)) {
    R <- cor(lead6p1, block6[, column])^2
    R <- t(R)
    blockresults6p1 <- cbind(blockresults6p1, R)
}
blockresults6p1 <- blockresults6p1[, -1]
blockresults6p1 <- t(blockresults6p1)
blockresults6p1 <- as.data.frame(blockresults6p1)
rownames(blockresults6p1) <- colnames(block6)
blockresults6p1$LDblock <- 6
colnames(blockresults6p1) <- c("R2", "LDblock")
for (column in 1:ncol(block6)) {
    R <- cor(lead6p2, block6[, column])^2
    R <- t(R)
    blockresults6p2 <- cbind(blockresults6p2, R)
}
blockresults6p2 <- blockresults6p2[, -1]
blockresults6p2 <- t(blockresults6p2)
blockresults6p2 <- as.data.frame(blockresults6p2)
rownames(blockresults6p2) <- colnames(block6)
blockresults6p2$LDblock <- 6
colnames(blockresults6p2) <- c("R2", "LDblock")
for (column in 1:ncol(block79)) {
    R <- cor(lead79p1, block79[, column])^2
    R <- t(R)
    blockresults79p1 <- cbind(blockresults79p1, R)
}
blockresults79p1 <- blockresults79p1[, -1]
blockresults79p1 <- t(blockresults79p1)
blockresults79p1 <- as.data.frame(blockresults79p1)
rownames(blockresults79p1) <- colnames(block79)
blockresults79p1$LDblock <- 79
colnames(blockresults79p1) <- c("R2", "LDblock")
for (column in 1:ncol(block79)) {
    R <- cor(lead79p2, block79[, column])^2
    R <- t(R)
    blockresults79p2 <- cbind(blockresults79p2, R)
}
blockresults79p2 <- blockresults79p2[, -1]
blockresults79p2 <- t(blockresults79p2)
blockresults79p2 <- as.data.frame(blockresults79p2)
rownames(blockresults79p2) <- colnames(block79)
blockresults79p2$LDblock <- 79
colnames(blockresults79p2) <- c("R2", "LDblock")

blockresults6 <- cbind(blockresults6p1[!is.na(blockresults6p1$R2),], blockresults6p2[!is.na(blockresults6p2$R2),])
blockresults79 <- cbind(blockresults79p1[!is.na(blockresults79p1$R2),], blockresults79p2[!is.na(blockresults79p2$R2),])
blockresults6 <- blockresults6[,-4]
blockresults79 <- blockresults79[,-4]
colnames(blockresults6) <- c("R21", "LDblock", "R22")
colnames(blockresults79) <- c("R21", "LDblock", "R22")
blockresults6$R2 <- pmax(blockresults6$R21, blockresults6$R22)
blockresults79$R2 <- pmax(blockresults79$R21, blockresults79$R22)
blockresults6 <- blockresults6[, c(2, 4)]
blockresults79 <- blockresults79[, c(2, 4)]

rm(block6, block79, lead6p1, lead6p2, lead79p1, lead79p2, blockresults6p1, blockresults6p2, blockresults79p1, blockresults79p2)


#Now, run the rest of the blocks through the correlation function
genoinput <- genoinput[genoinput$LDblock != 6 & genoinput$LDblock != 58 & genoinput$LDblock != 79,]
results <- data.frame(R2 = numeric(), LDblock = numeric())
for (block in unique(genoinput$LDblock)) {
    region <- genoinput[genoinput$LDblock == block,]
    lead <- region[region$lead == 1, 3:717]
    lead <- t(lead)
    nonlead <- region[, 3:717]
    nonlead <- t(nonlead)
    blockresults <- data.frame(2)
    for (column in 1:ncol(nonlead)) {
        R <- cor(lead, nonlead[, column])^2
        R <- t(R)
        blockresults <- cbind(blockresults, R)
    }
    blockresults <- blockresults[, -1]
    colnames(blockresults) <- colnames(nonlead)
    blockresults <- t(blockresults)
    blockresults <- as.data.frame(blockresults)
    blockresults$LDblock <- block
    colnames(blockresults) <- c("R2", "LDblock")
    results <- rbind(results, blockresults[!is.na(blockresults$R2), ])
}

#add the results from blocks 6 and 79
results <- rbind(results, blockresults6)
results <- rbind(results, blockresults79)
write.csv(results, "LDresults.csv", row.names = TRUE)

#Subset to only those which pass an R2 of 0.8 threshold
rsquarethreshold <- results[results$R2 >= 0.8,]
write.csv(rsquarethreshold, "R20.8.csv", row.names = TRUE)
#Tidy the workspace
rm(results, blockresults6, blockresults79, blockresults, genoinput, lead, nonlead, region, R, block, column)

#Load in the expression data for SN and filter for Europeans and those without Alzheimer's
SNexpression <- read.table("SNexpression.bed", sep = "\t", header = TRUE)
SNexpressioneur <- cbind(SNexpression[,1:4], SNexpression[, which(colnames(SNexpression) %in% phenotypes$SUBJID)])
phenosbrain <- phenotypes[phenotypes$SUBJID %in% colnames(SNexpressioneur),]
table(phenosbrain$MHALZDMT)
table(phenosbrain$MHPRKNSN)
table(phenosbrain$MHSMKSTS)
phenosbrain <- phenosbrain[phenosbrain$MHALZDMT == 0,]
gene_id <- SNexpressioneur$gene_id
SNexpressioneur <- cbind(gene_id, SNexpressioneur[, which(colnames(SNexpressioneur) %in% phenosbrain$SUBJID)])
rm(gene_id)
write.table(SNexpressioneur, "convertedSNexpression.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
SNexpressioneurcov <- read.table("convertedSNexpression.txt", sep = "\t", header = TRUE, row.names = "gene_id")
cisgenes <- read.csv("Ensemblgenes1.5mb.csv", header = TRUE)
rownames(SNexpressioneurcov) <- sub("\\..*", "", rownames(SNexpressioneurcov))
SNexpressioncis <- SNexpressioneurcov[rownames(SNexpressioneurcov) %in% cisgenes$name2,]
write.table(SNexpressioncis, "convertedSNexpressioncis.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

#Check to see if time since death is a confounder
class(SNexpressioneurcov)
SNexpressioneurcov <- as.matrix(SNexpressioneurcov)
phenosbrain <- subset(phenosbrain, select = c("SUBJID", "TRISCHD"))
phenosbrain$TRISCHD=as.numeric(phenosbrain$TRISCHD)
length(intersect(colnames(SNexpressioneurcov), phenosbrain$SUBJID))  
design=model.matrix(~TRISCHD, data = phenosbrain)

vfit <- lmFit(SNexpressioneurcov, design)
testfit <- lmFit(SNexpressioneurcov[1,], design)
names(testfit)
summary(testfit)
class(testfit)
vfite <- eBayes(vfit)
soptions(digits=3)
summary(vfit)
min(vfit$p.value)
sig <- vfit$p.value < 0.05
class(sig)
sig <- as.data.frame(sig)
table(sig$TRISCHD)
prop.table(table(sig$TRISCHD))



designsex <- model.matrix(~SEX + AGE, data = phenosbrain)
agesexfit <- lmFit(SNexpressioneurcov, designsex)
agesexfite <- eBayes(agesexfit)
names(agesexfit)
names(agesexfite)
testagesex <- lmFit(SNexpressioneurcov[1,], designsex)
testagesex$coefficients
names(testagesexe)
names(agesexfite)
testagesexe <- eBayes(testagesex)
testagesexe$var.p
sigage <- agesexfit$p.value < 0.05
#It is, so we need to adjust for it in the covariates
rm(vfit, Rownames, Matrix3, Matrix, design, sig)
phenosbrain <- subset(phenosbrain, select = c("SUBJID", "SEX", "AGE", "TRISCHD", "MHSMKSTS"))
covariates <- phenosbrain
covariates$MHSMKSTS <- ifelse(covariates$MHSMKSTS == "Yes", 1, 0)
covariates <- t(covariates)

covariates <- read.table("covariates.txt", sep = "\t", header = TRUE, row.names = 1)
sampleinfo <- read.csv("GTExsampleinfo.csv", header = TRUE)
sampleinfo <- sampleinfo[sampleinfo$SMTSD == "Brain - Substantia nigra",]
sampleinfo$SAMPID <- sub("-0011.*", "", sampleinfo$SAMPID)
sampleinfo$SAMPID <- sub("-", ".", sampleinfo$SAMPID)
covariates <- t(covariates)
sampleinfo <- sampleinfo[sampleinfo$SAMPID %in% rownames(covariates),]
sampleinfo <- sampleinfo[, c(1, 5)]
sampleinfo <- unique(sampleinfo)
sampleinfo <- sampleinfo[order(sampleinfo$SAMPID),]
covariates <- covariates[order(rownames(covariates)),]
covariates <- as.data.frame(covariates)
covariates$RIN <- sampleinfo$SMRIN
covariates <- t(covariates)
covariates$SUBJID <- rownames(covariates)
covariates <- covariates[, c(ncol(covariates), 1:(ncol(covariates)-1))]
write.table(covariates, "covariates.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
cova <- read.table("covariates.txt")
#Now, run the matrix eQTL package to find the eQTLs
#First, we need to make the files in the correct format
#convert the expression file to a format that can be used by the MatrixEQTL package
common_cols <- intersect(colnames(genoinfo), colnames(SNexpressioneur))
convertedgenotypeSN <- genoinfo[, "ID"]
convertedgenotypeSN <- cbind(convertedgenotypeSN, genoinfo[, common_cols])
write.table(convertedgenotypeSN, "PDlociSNb38.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)





#Set required configuration for the MatrixEQTL package
blockcoords <- read.table("blockcoordsb38.bed", header = TRUE, sep = "\t")
blockcoords$start <- blockcoords$chromStart - 1500000
blockcoords$end <- blockcoords$chromEnd + 1500000
blockcoords <- blockcoords[, c(1, 4, 5)]
write.table(blockcoords, "b381.5mb.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
cisgenesgtex <- read.csv("1.5mbcisgenesb38.csv", header = TRUE)
cisgenes <- read.csv("Ensemblgenes1.5mb.csv")
expression <- read.table("convertedSNexpression.txt", header = TRUE, sep = "\t")
length(intersect(cisgenesgtex$geneId, cisgenes$name2))
length(intersect(cisgenesgtex$geneId, expression$gene_id))
cisexpression <- expression[expression$gene_id %in% cisgenesgtex$geneId,]
cisgenesgtexexpressed <- cisgenesgtex[cisgenesgtex$geneId %in% cisexpression$gene_id,]
cisgenesgtexexpressed <- unique(cisgenesgtexexpressed)
blockcoords$LDblock <- c(1:88)
colnames(cisgenesgtexexpressed) <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand", "geneId", "geneType", "expCount", "expScores")
blockgenes <- left_join(blockcoords, cisgenesgtexexpressed, by = "chrom")
blockgenes <- blockgenes[blockgenes$chromStart >= blockgenes$start & blockgenes$chromStart <= blockgenes$end,]
blockgenes <- unique(blockgenes)
length(unique(blockgenes$geneId))
length(unique(cisgenesgtexexpressed$geneId))
table(blockgenes$LDblock)
blocksumm <- blockgenes %>%
    group_by(LDblock) %>%
    summarise(num_genes = n_distinct(geneId),
    )
write.csv(blocksumm, "blockgenesummary.csv", row.names = FALSE)

#Now, create the files for the eQTL analysis
for (block in c(2:57, 59:88)) {
    genes <- blockgenes[blockgenes$LDblock == block,]
    blockexp <- expression[expression$gene_id %in% genes$geneId,]
    write.table(blockexp, paste("block", block, "exp.txt", sep = ""), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    r2snps <- read.csv("R20.8.csv", row.names = 1)
    r2snps <- r2snps[r2snps$LDblock == block,]
    LDresults <- read.csv("LDresults.csv", row.names = 1)
    LDresults <- LDresults[LDresults$LDblock == block,]
    geno <- read.table("PDlociSNb38.txt", header = TRUE, sep = "\t")
    geno <- geno[geno$convertedgenotypeSN %in% rownames(r2snps),]
    write.table(geno, paste("block", block, "genoall.txt", sep = ""), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    r2snps <- r2snps[r2snps$R2 == 1,]
    genolead <- geno[geno$convertedgenotypeSN %in% rownames(r2snps),]
    write.table(genolead, paste("block", block, "genolead.txt", sep = ""), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    allsnps <- read.csv("LDresults.csv", row.names = 1)
    allsnps <- allsnps[allsnps$LDblock == block,]
    genonoLD <- read.table("PDlociSNb38.txt", header = TRUE, sep = "\t")
    genonoLD <- genonoLD[genonoLD$convertedgenotypeSN %in% rownames(allsnps),]
    write.table(genonoLD, paste("block", block, "genonoLD.txt", sep = ""), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}





covariates_file_name = "covariates.txt"
#Load in the covariates
cvrt = SlicedData$new();
cvrt$fileDelimiter = '\t'; # the TAB character
cvrt$fileOmitCharacters = 'NA'; # denote missing values;
cvrt$fileSkipRows = 1; # one row of column labels
cvrt$fileSkipColumns = 1; # one column of row labels
cvrt$fileSliceSize = snps$nCols()+1; # read file in one piece
if(length(covariates_file_name)>0) {
    cvrt$LoadFile(covariates_file_name);
}

setwd("C:/Users/Thomas Brightwell/OneDrive - Imperial College London/PD_MSc_Project2024/Working materials/Visual Studio/eQTL/data")
for (block in c(2:57, 59:88)) {
    base.dir = find.package("MatrixEQTL")
    useModel = modelLINEAR
    SNP_file_name = paste("block", block, "genoall.txt", sep = "")
    expression_file_name = paste("block", block, "exp.txt", sep = "")
    output_file_name = paste("block", block, "eqtlallv2.txt", sep = "")
    pvOutputThreshold = 0.05
    errorCovariance = numeric()


    #Load in the genotype data
    snps = SlicedData$new()
    snps$fileDelimiter = "\t"
    snps$fileOmitCharacters = "NA"
    snps$fileSkipRows = 1
    snps$fileSkipColumns = 1
    snps$fileSliceSize = 2000
    snps$LoadFile( SNP_file_name )

    #Load in the expression data
    gene = SlicedData$new();
    gene$fileDelimiter = '\t'; # the TAB character
    gene$fileOmitCharacters = 'NA'; # denote missing values;
    gene$fileSkipRows = 1; # one row of column labels
    gene$fileSkipColumns = 1; # one column of row labels
    gene$fileSliceSize = 10000; # read file in pieces of 10,000 rows
    gene$LoadFile(expression_file_name);


    #Run the eQTL analysis
    SNeQTL = Matrix_eQTL_engine(
        snps = snps,
        gene = gene,
        cvrt = cvrt,
        output_file_name = output_file_name,
        pvOutputThreshold = pvOutputThreshold,
        useModel = useModel,
        errorCovariance = errorCovariance,
        verbose = TRUE,
        pvalue.hist = TRUE,
        min.pv.by.genesnp = FALSE,
        noFDRsaveMemory = FALSE
    )

    base.dir = find.package("MatrixEQTL")
    useModel = modelLINEAR
    SNP_file_name = paste("block", block, "genolead.txt", sep = "")
    expression_file_name = paste("block", block, "exp.txt", sep = "")
    output_file_name = paste("block", block, "eqtlleadv2.txt", sep = "")
    pvOutputThreshold = 0.05
    errorCovariance = numeric()


    #Load in the genotype data
    snps = SlicedData$new()
    snps$fileDelimiter = "\t"
    snps$fileOmitCharacters = "NA"
    snps$fileSkipRows = 1
    snps$fileSkipColumns = 1
    snps$fileSliceSize = 2000
    snps$LoadFile( SNP_file_name )

    #Load in the expression data
    gene = SlicedData$new();
    gene$fileDelimiter = '\t'; # the TAB character
    gene$fileOmitCharacters = 'NA'; # denote missing values;
    gene$fileSkipRows = 1; # one row of column labels
    gene$fileSkipColumns = 1; # one column of row labels
    gene$fileSliceSize = 10000; # read file in pieces of 10,000 rows
    gene$LoadFile(expression_file_name);


    #Run the eQTL analysis
    SNeQTL = Matrix_eQTL_engine(
        snps = snps,
        gene = gene,
        cvrt = cvrt,
        output_file_name = output_file_name,
        pvOutputThreshold = pvOutputThreshold,
        useModel = useModel,
        errorCovariance = errorCovariance,
        verbose = TRUE,
        pvalue.hist = TRUE,
        min.pv.by.genesnp = FALSE,
        noFDRsaveMemory = FALSE
    )
     base.dir = find.package("MatrixEQTL")
    useModel = modelLINEAR
    SNP_file_name = paste("block", block, "genonoLD.txt", sep = "")
    expression_file_name = paste("block", block, "exp.txt", sep = "")
    output_file_name = paste("block", block, "eqtlnoLDv2.txt", sep = "")
    pvOutputThreshold = 0.05
    errorCovariance = numeric()


    #Load in the genotype data
    snps = SlicedData$new()
    snps$fileDelimiter = "\t"
    snps$fileOmitCharacters = "NA"
    snps$fileSkipRows = 1
    snps$fileSkipColumns = 1
    snps$fileSliceSize = 2000
    snps$LoadFile( SNP_file_name )

    #Load in the expression data
    gene = SlicedData$new();
    gene$fileDelimiter = '\t'; # the TAB character
    gene$fileOmitCharacters = 'NA'; # denote missing values;
    gene$fileSkipRows = 1; # one row of column labels
    gene$fileSkipColumns = 1; # one column of row labels
    gene$fileSliceSize = 10000; # read file in pieces of 10,000 rows
    gene$LoadFile(expression_file_name);


    #Run the eQTL analysis
    SNeQTL = Matrix_eQTL_engine(
        snps = snps,
        gene = gene,
        cvrt = cvrt,
        output_file_name = output_file_name,
        pvOutputThreshold = pvOutputThreshold,
        useModel = useModel,
        errorCovariance = errorCovariance,
        verbose = TRUE,
        pvalue.hist = TRUE,
        min.pv.by.genesnp = FALSE,
        noFDRsaveMemory = FALSE
    )
}


#bind the results together into one file each
eqtlall <- data.frame()
eqtllead <- data.frame()
eqtlnoLD <- data.frame()
for (block in c(2:57, 59:88)) {
    res <- read.table(paste("block", block, "eqtlallv2.txt", sep = ""), header = TRUE, sep = "\t")
    eqtlall <- rbind(eqtlall, res)
    reslead <- read.table(paste("block", block, "eqtlleadv2.txt", sep = ""), header = TRUE, sep = "\t")
    eqtllead <- rbind(eqtllead, reslead)
    resnoLD <- read.table(paste("block", block, "eqtlnoLDv2.txt", sep = ""), header = TRUE, sep = "\t")
    eqtlnoLD <- rbind(eqtlnoLD, resnoLD)
}


write.table(eqtllead, "eqtllead.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

write.table(eqtlnoLD, "eqtlnoLD.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

blockgenes <- blockgenes[blockgenes$LDblock !=1,]

length(unique(eqtl$gene))
eqtlall$LDblock <- blockgenes$LDblock[match(eqtlall$gene, blockgenes$geneId)]
table(eqtlall$LDblock)
length(unique(eqtlall$LDblock))
eqtlsumm <- eqtlall %>%
    group_by(LDblock) %>%
    summarise(num_genes = n_distinct(gene))
write.csv(eqtlsumm, "eqtlallsummary.csv", row.names = FALSE)
write.table(eqtlall, "eqtlallv2.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

eqtllead$LDblock <- blockgenes$LDblock[match(eqtllead$gene, blockgenes$geneId)]
table(eqtllead$LDblock)
length(unique(eqtllead$LDblock))
eqtlleadsumm <- eqtllead %>%
    group_by(LDblock) %>%
    summarise(num_genes = n_distinct(gene))
write.csv(eqtlleadsumm, "eqtlleadsummary.csv", row.names = FALSE)
eqtlnoLD$LDblock <- blockgenes$LDblock[match(eqtlnoLD$gene, blockgenes$geneId)]
table(eqtlnoLD$LDblock)
length(unique(eqtlnoLD$LDblock))
eqtlnoLDsumm <- eqtlnoLD %>%
    group_by(LDblock) %>%
    summarise(num_genes = n_distinct(gene))
write.csv(eqtlnoLDsumm, "eqtlnoLDsummary.csv", row.names = FALSE)


length(unique(eqtl$gene[eqtl$FDR <= 0.05])) 
length(unique(eqtllead$gene[eqtllead$FDR <= 0.05]))
length(unique(eqtlnoLD$gene[eqtlnoLD$FDR <= 0.05]))
table(eqtlnoLD$LDblock[eqtlnoLD$FDR <= 0.05])
#add a column to indicate if the gene is mitochondrial or not
mitocarta <- read.csv("mitocarta.csv", header = TRUE)
genes$mitochondrial <- ifelse(genes$gene %in% mitocarta$EnsemblGeneID_mapping_version_20200130, 1, 0)
table(genes$mitochondrial)
write.table(genes, "eQTLgenesmito.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

genes <- read.table("eQTLgenesmito.txt", sep = "\t", header = TRUE)
test <- read.GMT("testgmt.gmt")
head(test)
allcisgenes <- genes$gene
allcisgenes <- as.list(allcisgenes)
class(testgenes)
is.GMT(testgenes)
write.GMT(testgenes, "testgenes.gmt")
writeGMT <- function #Create a gmt (gene matrix transposed) file
### Createss a gmt (gene matrix transposed) file such as those
### provided by mSigDB or geneSigDB, from an R list object.
### Function by Levi Waldron.
(object,
### R list object that will be converted to GMT file.  Each element
### should contain a vector of gene names, and the names of the
### elements will used for the gene set names
 fname
### Output file name for .gmt file
 ){
  if (class(object) != "list") stop("object should be of class 'list'")
  if(file.exists(fname)) unlink(fname)
  for (iElement in 1:length(object)){
    write.table(t(c(make.names(rep(names(object)[iElement],2)),object[[iElement]])),
                sep="\t",quote=FALSE,
                file=fname,append=TRUE,col.names=FALSE,row.names=FALSE)
  }
### Called for the effect of writing a .gmt file
}
###https://github.com/lwaldron/LeviRmisc/blob/master/R/writeGMT.R
writeGMT(allcisgenes, "cisgenes.gmt")
mitocisgenes <- genes[genes$mitochondrial == 1,]
mitocisgenes <- mitocisgenes$gene
mitocisgenes <- as.list(mitocisgenes)
writeGMT(mitocisgenes, "mitocisgenes.gmt")

zscores <- read.csv("Zscores.csv", header = TRUE)
colnames(zscores) <- c("gene", "unknown", "zscore", "pvalue", "pvalue2", "adjpvalue", "adjpvalue2")
zscores <- zscores[zscores$gene %in% mitocarta$Symbol,]
write.csv(zscores, "mitozscores.csv", row.names = FALSE)





eqtl_summary <- eqtl %>%
    group_by(LDblock) %>%
    summarise(num_eqtls = n(), unique_genes = n_distinct(gene), unique_snps = n_distinct(SNP))











convert.vcf("PDlocib38.recode.vcf", which = "GT", map = TRUE, snp.pos = FALSE, genotype_file_name = "PDlocib38.txt")
genotypeb38short <- read.table("PDlocib38.txt", sep = "\t", header = TRUE)
genotypeb38short <- inner_join(genotypeb38short, lookuptable, join_by("snpid" == "variant_id"))
ID <- subset(genotypeb38short, select = c("snpid"))
ID$lead <- ifelse(ID$snpid %in% leadSNPs$variant_id, 1, 0)
present <- anti_join(leadSNPs, ID, join_by("variant_id" == "snpid"))





test2 <- genotypeb38short[genotypeb38short$variant_id == "rs114",]







samples <- genotypeb38short[, 2:839]
info <- genotypeb38short[, 840:846]
genotypeb38short <- cbind(ID, info, samples)
rm(ID, samples, info)
write.table(genotypeb38short, "annotated genotype for PDloci.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


#Conduct regression to calculate LD statistics
genotypeb38short <- read.table("annotated genotype for PDloci.txt", sep = "\t", header = TRUE)





#### Specifically for SNCA, rerun an analysis.
SNCAregion <- read.table("chr4.recode.vcf", sep = "\t", header = TRUE)
SNCAregion <- SNCAregion[SNCAregion$CHROM == "chr4" & SNCAregion$POS >= 89685444 & SNCAregion$POS <= 89856062,]
SNCAregion <- SNCAregion[, colnames(SNCAregion) %in% colnames(genoinfo)]
for (col in 10:ncol(SNCAregion)) {
    SNCAregion[, col] <- sub(":.*", "", SNCAregion[, col])
    SNCAregion[, col] <- gsub("0/0", "0", SNCAregion[, col])
    SNCAregion[, col] <- gsub("0/1", "1", SNCAregion[, col])
    SNCAregion[, col] <- gsub("1/1", "2", SNCAregion[, col])
    SNCAregion[, col] <- gsub("1/0", "1", SNCAregion[, col])
    SNCAregion[, col] <- gsub("\\./\\.", "3", SNCAregion[, col])
    SNCAregion[, col] <- as.numeric(SNCAregion[, col])
    SNCAregion[, col][SNCAregion[, col] == 3] <- NA
}
rownames(SNCAregion) <- SNCAregion$ID

blocks28and29 <- SNCAregion[, 10:ncol(SNCAregion)]
lead28 <- blocks28and29[which(rownames(blocks28and29) == "chr4_89704960_G_A_b38"),]
lead29 <- genoinput[which(rownames(genoinput) == "chr4_89715557_T_C_b38"), 3:ncol(genoinput)]
lead28 <- t(lead28)
lead29 <- t(lead29)
blocks28and29 <- t(blocks28and29)
blockresults28 <- data.frame(2)
blockresults29 <- data.frame(2)
for (column in 1:ncol(blocks28and29)) {
    R <- cor(lead28, blocks28and29[, column], use = "complete.obs")^2
    R <- t(R)
    blockresults28 <- cbind(blockresults28, R)
}
blockresults28 <- blockresults28[, -1]
blockresults28 <- t(blockresults28)
blockresults28 <- as.data.frame(blockresults28)
rownames(blockresults28) <- colnames(blocks28and29)
blockresults28$LDblock <- 28
colnames(blockresults28) <- c("R2", "LDblock")
for (column in 1:ncol(blocks28and29)) {
    R <- cor(lead29, blocks28and29[, column], use = "complete.obs")^2
    R <- t(R)
    blockresults29 <- cbind(blockresults29, R)
}
blockresults29 <- blockresults29[, -1]
blockresults29 <- t(blockresults29)
blockresults29 <- as.data.frame(blockresults29)
rownames(blockresults29) <- colnames(blocks28and29)
blockresults29$LDblock <- 29
colnames(blockresults29) <- c("R2", "LDblock")
blockresults28and29 <- cbind(blockresults28[!is.na(blockresults28$R2),], blockresults29[!is.na(blockresults29$R2),])
colnames(blockresults28and29) <- c("R21", "LDblock", "R22", "LDblock2")
blockresults28and29$R2 <- pmax(blockresults28and29$R21, blockresults28and29$R22)
blockresults28and29 <- blockresults28and29[, c(2, 5)]
blockresults28and29$ID <- rownames(blockresults28and29)

write.table(blockresults28and29, "LDresults28and29.csv", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
nold2829 <- read.table("LDresults28and29.csv", sep = "\t", header = TRUE)

r2threshold <- blockresults28and29[blockresults28and29$R2 >= 0.8,]
blocks28and29 <- SNCAregion[, 10:ncol(SNCAregion)]
blocks28and29 <- blocks28and29[rownames(blocks28and29) %in% r2threshold$ID,]
blocks28and29$ID <- rownames(blocks28and29)
blocks28and29 <- blocks28and29[, c(ncol(blocks28and29), 1:(ncol(blocks28and29)-1))]
block29 <- read.table("block29genoall.txt", sep = "\t", header = TRUE)
colnames(block29)[1] <- "ID"
blocks28and29 <- blocks28and29[,colnames(blocks28and29) %in% colnames(block29)]
write.table(blocks28and29, "block28and29genoall.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


blocks28and29nold <- blocks28and29[rownames(blocks28and29) %in% nold2829$ID,]
blocks28and29nold$ID <- rownames(blocks28and29nold)
blocks28and29nold <- blocks28and29nold[, c(ncol(blocks28and29nold), 1:(ncol(blocks28and29nold)-1))]
blocks28and29nold <- blocks28and29nold[,colnames(blocks28and29nold) %in% colnames(block29)]
blocks28and29nold <- blocks28and29nold[, colnames(blocks28and29nold) %in% colnames(block29)]
write.table(blocks28and29nold, "block28and29genonoLD.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)



base.dir = find.package("MatrixEQTL")
useModel = modelLINEAR
SNP_file_name = "block28and29genonoLD.txt"
expression_file_name = "block28exp.txt"
output_file_name = "SNCAeqtlnoLD.txt"
covariates_file_name = "covariates.txt"
pvOutputThreshold = 0.05
errorCovariance = numeric()


#Load in the genotype data
snps = SlicedData$new()
snps$fileDelimiter = "\t"
snps$fileOmitCharacters = "NA"
snps$fileSkipRows = 1
snps$fileSkipColumns = 1
snps$fileSliceSize = 2000
snps$LoadFile( SNP_file_name )

#Load in the expression data
gene = SlicedData$new();
gene$fileDelimiter = '\t'; # the TAB character
gene$fileOmitCharacters = 'NA'; # denote missing values;
gene$fileSkipRows = 1; # one row of column labels
gene$fileSkipColumns = 1; # one column of row labels
gene$fileSliceSize = 10000; # read file in pieces of 10,000 rows
gene$LoadFile(expression_file_name);



#Load in the covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = '\t'; # the TAB character
cvrt$fileOmitCharacters = 'NA'; # denote missing values;
cvrt$fileSkipRows = 1; # one row of column labels
cvrt$fileSkipColumns = 1; # one column of row labels
cvrt$fileSliceSize = snps$nCols()+1; # read file in one piece
if(length(covariates_file_name)>0) {
    cvrt$LoadFile(covariates_file_name);
}

#Run the eQTL analysis
SNeQTL = Matrix_eQTL_engine(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name = output_file_name,
    pvOutputThreshold = pvOutputThreshold,
    useModel = useModel,
    errorCovariance = errorCovariance,
    verbose = TRUE,
    pvalue.hist = TRUE,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE
)

eqtlSNCA <- read.table("SNCAeqtl.txt", sep = "\t", header = TRUE)
eqtlSNCA <- eqtlSNCA[order(eqtlSNCA$p.value),]
eqtlSNCA <- eqtlSNCA[!duplicated(eqtlSNCA$gene),]


eqtlSNCAnoLD <- read.table("SNCAeqtlnoLD.txt", sep = "\t", header = TRUE)
eqtlSNCAnoLD <- eqtlSNCAnoLD[order(eqtlSNCAnoLD$p.value),]
eqtlSNCAnoLD <- eqtlSNCAnoLD[eqtlSNCAnoLD$gene == "ENSG00000145335.15",]
eqtlSNCAnoLD$R2 <- nold2829$R2[match(eqtlSNCAnoLD$SNP, nold2829$ID)]
