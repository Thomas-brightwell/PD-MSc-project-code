setwd("C:/Users/Thomas Brightwell/OneDrive - Imperial College London/PD_MSc_Project2024/Working materials/Visual Studio/BRAINEAC")
first20 <- read.table("BRAINEACsnps1-20.tsv", sep = "\t", header = TRUE)
second20 <- read.table("BRAINEACsnps21-40.tsv", sep = "\t", header = TRUE)
third20 <- read.table("BRAINEACsnps41-60.tsv", sep = "\t", header = TRUE)
fourth20 <- read.table("BRAINEACsnps61-80.tsv", sep = "\t", header = TRUE)
last10 <- read.table("BRAINEACsnps81-90.tsv", sep = "\t", header = TRUE)
allresults <- rbind(first20, second20, third20, fourth20, last10)
rm(first20, second20, third20, fourth20, last10)
leadsnps <- subset(allresults, select = c(rsid, chr))
leadsnps <- unique(leadsnps)

install.packages("tidyverse")
library(tidyverse)
leadsnps <- arrange(leadsnps)
table2 <-read.table("TableS2_v2.txt", header = TRUE, sep = "\t")
table2 <- select(table2, snp, chr)
colnames(table2) <- c("rsid", "chr")
missingsnps <- anti_join(table2, leadsnps, by = "rsid")
write.table(missingsnps, "missingsnps.tsv", sep = "\t", row.names = FALSE)
#import the data from a neighbouring SNP for the missing lead SNPs





write.table(allresults, "allresults.tsv", sep = "\t", row.names = FALSE)
View(allresults)
signif <- subset(allresults, allresults$SNIG <= 0.05)
write.table(signif, "significant_transcripts_SN.tsv", sep = "\t", row.names = FALSE)
genes <- unique(signif$geneSymbol)
write.table(genes, "signif_genes_SN.tsv", sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
