#This script compares the overlap of eQTL disease cis-genes with the DGE analysis results
setwd("C:/Users/Thomas Brightwell/OneDrive - Imperial College London/PD_MSc_Project2024/DGE & GSEA")
library(tidyverse)
BiocManager::install("fgsea")
library(fgsea)
install.packages("ggExtra")
library(ggExtra)
library(biomaRt)
library(data.table)
if (!require(devtools)) install.packages("devtools")
devtools::install_github("nicolash2/ggvenn")
library(ggvenn)
install.packages("UpSetR")
library(UpSetR)
install.packages("ComplexUpset")
library(ComplexUpset)
#Load the GTEx eQTL data
GTEx <- read.table("eqtlallv2.txt", sep = "\t", header = TRUE)
GTExsymbols <- read.table("gtexgenesymbols.txt", sep = ",", header = TRUE)
GTEx$gene <- sub("\\..*", "", GTEx$gene)
GTEx <- merge(GTEx, GTExsymbols, by.x = "gene", by.y = "Gene.stable.ID", all.x = TRUE)
GTEx <- GTEx[order(GTEx$p.value),]
GTEx <- GTEx[!duplicated(GTEx$gene),]
GTEx <- GTEx[, c(2, 8, 5, 7)]
colnames(GTEx) <- c("SNP", "gene", "p.value", "LDblock")
GTEx <- GTEx[GTEx$gene != "",]
rm(GTExsymbols)

#Load the BRAINEAC eQTL data
BRAINEAC <- read.table("BRAINEACeqtlPMI.txt", sep = "\t", header = TRUE)
BRAINEAC <- BRAINEAC[order(BRAINEAC$p.value),]
BRAINEAC <- BRAINEAC[!duplicated(BRAINEAC$gene),]
BRAINEAC <- BRAINEAC[, c(1, 2, 5, 7)]
eQTLboth <- read.table("overlapPMI.txt", sep = "\t", header = TRUE)

overlap <- GTEx[GTEx$gene %in% BRAINEAC$gene,]

#Combine the results into one list of genes
combined <- rbind(GTEx, BRAINEAC)




combined <- combined[complete.cases(combined),]
combined <- combined[!duplicated(combined$gene),]
write.table(combined, "combinedeqtls.txt", sep = "\t", row.names = FALSE, quote = FALSE)
combined <- read.table("combinedeqtls.txt", sep = "\t", header = TRUE)


#Load the DGE analysis results for GEO
GEO <- read.csv("GEOzscores.csv")
GEOcisgenes <- GEO[GEO$X %in% combined$gene,]
sigGEOcis <-GEOcisgenes[GEOcisgenes$p_values2 <= 0.05,]
#Load the DGE analysis results for PPMI
PPMI <- read.csv("PPMIzscores.csv")
PPMIcisgenes <- PPMI[PPMI$X %in% combined$gene,]
sigPPMIcis <-PPMIcisgenes[PPMIcisgenes$P.Value <= 0.05,]

eitheror <- c(sigGEOcis$X, sigPPMIcis$X)
eitheror <- unique(eitheror)
sigboth <- intersect(sigGEOcis$X, sigPPMIcis$X)

botheQTL <- read.table("overlapPMI.txt", sep = "\t", header = TRUE)
GEOsignif <- GEO[GEO$p_values2 <= 0.05,]
PPMIsignif <- PPMI[PPMI$P.Value <= 0.05,]

bestdata <- intersect(eitheror, GTEx$gene)




write.table(eitheror, "eitheror.txt", sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)

minilist <- Reduce(intersect, list(GEO$X, PPMI$X, counts$X, GTExexprnames$hgnc_symbol))


minilist <- Reduce(intersect, list(GEO$X[GEO$p_values2 <= 0.05], BRAINEACnoLD$gene, GTEx$gene))

test6 <- intersect(GEO$X[GEO$p_values2 <= 0.05], PPMI$X[PPMI$P.Value <= 0.05])

sharedall <- Reduce(intersect, list(GEOsignif$X, combined$gene))





GNE <- intersect(GEO$X[GEO$p_values2 <= 0.05], GTEx$Gene.name)
GNP <- intersect(PPMI$X[PPMI$P.Value <= 0.05], GTEx$Gene.name)
BNE <- intersect(GEO$X[GEO$p_values2 <= 0.05], BRAINEAC$gene)
BNP <- intersect(PPMI$X[PPMI$P.Value <= 0.05], BRAINEAC$gene)

unionall <- c(GNE, GNP, BNE, BNP)
unionall <- unique(unionall)
length(unionall)


gtexeqtlsymbols <- GTEx[GTEx$Gene.name != "",]
gtexeqtlsymbols <- GTEx[GTEx$Gene.name]#
write.table(GTEx, "GTExeqtls.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(BRAINEAC, "BRAINEACeqtls.txt", sep = "\t", row.names = FALSE, quote = FALSE)




#assemble a list of cisgenes tested

regionsb37 <- read.table("1.5mbcoordsb37.bed", sep = "\t", header = FALSE)
regionsb37$V1 <- paste0("chr", regionsb37$V1)
colnames(regionsb37) <- c("Chr", "Start37", "End37")
regionsb37$block <- c(1:88)
regionsb38 <- read.table("b381.5mb.bed", sep = "\t", header = FALSE)
colnames(regionsb38) <- c("Chr", "Start38", "End38")
regionsb38$block <- c(1:88)
regions <- merge(regionsb37, regionsb38, by = c("Chr", "block"))
rm(regionsb37, regionsb38)


braineacexp <- read.csv("CountTable_69_SN.csv")
gtexexpr <- read.table("SNexpression.bed", sep = "\t", header = TRUE)
gtexexpr$gene_id <- sub("\\..*", "", gtexexpr$gene_id)

braineacexp <- braineacexp[,1]
braineacexp <- as.data.frame(braineacexp)
colnames(braineacexp) <- c("Name")
ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")


braineacexpchrs <- getBM(attributes = c("chromosome_name"), filters = "external_gene_name", values = braineacexp$Name, mart = ensembl)

GTExexprnames <- getBM(attributes = c("hgnc_symbol"), filters = "ensembl_gene_id", values = gtexexpr$gene_id, mart = ensembl)

##################
#MITO TESTING

GTExexprnames$impi <- GTExexprnames$hgnc_symbol %in% impionly$Symbol
table(GTExexprnames$impi)

mitocarta <- read.csv("mitocarta.csv", header = TRUE)
mitocarta <- mitocarta[mitocarta$Symbol != "",]

impi <- read.csv("impi.csv", header = TRUE)
impionly <- impi[impi$IMPI.Class != "Unclassified",]
impionly <- impionly[impionly$IMPI.Class != "Probable contaminant",]
impionly <- impionly[impionly$IMPI.Class != "Not reported mitochondrial",]
impionly <- impionly[impionly$IMPI.Class != "Dubious gene",]
impionly <- impionly[impionly$IMPI.Class != "",]

cisgenes <- read.csv("1.5mbcisgenesb38.csv", header = TRUE)
cisgenes <- cisgenes[cisgenes$name %in% c(GTExexprnames$hgnc_symbol, braineacexp$X),]
cisgenes <- cisgenes[!duplicated(cisgenes$name),]
cisgenes$eQTL <- cisgenes$name %in% combined$gene
table(cisgenes$eQTL)

GOmito <- read.csv("mito.csv", header = FALSE)
length(unique(GOmito$V3))
GOmito <- GOmito[!duplicated(GOmito$V3),]
cisgenes <- cisgenes[!duplicated(cisgenes$name),]
eitheror <- as.data.frame(eitheror)
eitheror$mito <- eitheror$eitheror %in% mitocarta$Symbol
eitheror$mitoGO <- eitheror$eitheror %in% GOmito$V3
eitheror$impi <- eitheror$eitheror %in% impionly$Symbol
table(eitheror$mito)
table(eitheror$mitoGO)
table(eitheror$impi)

#####################

#making fisher's exact test tables


combined$mito <- combined$gene %in% mitocarta$Symbol
combined$mitoGO <- combined$gene %in% GOmito$V3
combined$impi <- combined$gene %in% impionly$Symbol
combined$cis <- combined$gene %in% cisgenes$name
cisgenes$mito <- cisgenes$name %in% mitocarta$Symbol
cisgenes$mitoGO <- cisgenes$name %in% GOmito$V3
cisgenes$impi <- cisgenes$name %in% impionly$Symbol
table(cisgenes$mito)
table(cisgenes$mitoGO)
table(cisgenes$impi)
table(combined$mito)
table(combined$mitoGO)
table(combined$impi)
cisgenesnoeqtlnomito <- length(cisgenes$name[cisgenes$eQTL == FALSE & cisgenes$impi == FALSE])
cisgenesnoeqtlmito <- length(cisgenes$name[cisgenes$eQTL == FALSE & cisgenes$impi == TRUE])
cigeneseqtlmito <- length(cisgenes$name[cisgenes$eQTL == TRUE & cisgenes$impi == TRUE])
cisgeneseqtlnomito <- length(cisgenes$name[cisgenes$eQTL == TRUE & cisgenes$impi == FALSE])

fisherdata <- data.frame(
  "No eQTL" = c(cisgenesnoeqtlnomito, cisgenesnoeqtlmito),
  "eQTL" = c(cisgeneseqtlnomito, cigeneseqtlmito),
  row.names = c("No mito", "Mito"),
  stringsAsFactors = FALSE
)

mosaicplot(fisherdata, color = c("red", "blue"), main = "Mosaic plot for eQTL and mitochondrial genes", xlab = "Mitochondrial gene", ylab = "eQTL", las = 1)
chisq.test(fisherdata)$expected
fisher <- fisher.test(fisherdata)
fisher

combined$DGE <- combined$gene %in% eitheror$eitheror
eqtldgemito <- length(combined$gene[combined$DGE == TRUE & combined$impi == TRUE])
eqtldgenomito <- length(combined$gene[combined$DGE == TRUE & combined$impi == FALSE])
eqtlnodgemito <- length(combined$gene[combined$DGE == FALSE & combined$impi == TRUE])
eqtlnodgenomito <- length(combined$gene[combined$DGE == FALSE & combined$impi == FALSE])

fisherdata2 <- data.frame(
  "No DGE" = c(eqtlnodgemito, eqtlnodgenomito),
  "DGE" = c(eqtldgemito, eqtldgenomito),
  row.names = c("Mito", "no mito"),
  stringsAsFactors = FALSE
)

mosaicplot(fisherdata2, color = c("red", "blue"), main = "Mosaic plot for DGE and mitochondrial genes", xlab = "Mitochondrial gene", ylab = "DGE", las = 1)

fisher2 <- fisher.test(fisherdata2)
fisher2




fisherdata3 <- data.frame(
  "DEG" = c(978, 7767),
  "DEG" = c(791, 9193),
  row.names = c("Mito", "no mito"),
  stringsAsFactors = FALSE
)

fisher3 <- fisher.test(fisherdata3)
fisher3




###########################
#GSEA 






#Run fgsea on the combined list of genes
mito43 <- gmtPathways("mito43.gmt")
GEOscore <- GEO[, c(3, 1)]
GEOscore <- setNames(GEOscore$Z, GEOscore$X)
GEOgsea <- fgsea(
  mito43,
  GEOscore,
  minSize = 5,
  maxSize = 500,
)
write.table(GEOscore, "GEOscore.txt", sep = "\t", row.names = TRUE, quote = FALSE)
GEOscore <- read.table("GEOscore.txt", sep = "\t", header = FALSE)
GEOscore <- setNames(GEOscore$V2, GEOscore$V1)



cluster_counts.matrix <- cluster_counts.matrix[rowSums(cluster_counts.matrix >= 3) >= 0.2 * ncol(cluster_counts.matrix), ]
PPMIscore <- PPMIscore[PPMIscore$X %in% rownames(cluster_counts.matrix), ]

PPMIscore <- PPMI[, c(2, 1)]
PPMIscore <- setNames(PPMIscore$t, PPMIscore$X)
PPMIgsea <- fgsea(
  mito43,
  PPMIscore,
  minSize = 5,
  maxSize = 500,
)
write.table(PPMIscore, "PPMIscore.txt", sep = "\t", row.names = TRUE, quote = FALSE)
PPMIscore <- read.table("PPMIscore.txt", sep = "\t", header = FALSE)
PPMIscore <- setNames(PPMIscore$V2, PPMIscore$V1)


topPathwaysUpGEO <- GEOgsea[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDownGEO <- GEOgsea[ES < 0][head(order(pval), n=20), pathway]
topPathwayGEO <- c(topPathwaysUpGEO, rev(topPathwaysDownGEO))
ordered <- GEOgsea[order(NES),]

mito43[["REACTOME_RESPIRATORY_ELECTRON_TRANSPORT_ATP_SYNTHESIS_BY_CHEMIOSMOTIC_COUPLING_AND_HEAT_PRODUCTION_BY_UNCOUPLING_PROTEINS"]] <- "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT_ATP_SYNTHESIS"
GEOgsea[38,1] <- "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT_ATP_SYNTHESIS"
png("FullGEOGSEAplot.png", width = 800, height = 900)
plotGseaTable(mito43, GEOscore, ordered, pathwayLabelStyle = list(size=8), gseaParam=0.5, colwidths = c(6.2, 3, 0.8, 1.2, 1.2))
dev.off()
fwrite(GEOgsea, "GEOgsea.csv", sep=",", sep2=c("", " ", ""))
fwrite(PPMIgsea, "PPMIgsea.csv", sep=",", sep2=c("", " ", ""))


GEOgsea <- read.csv("GEOgsea.csv", header = TRUE)
PPMIgsea <- read.csv("PPMIgsea.csv", header = TRUE)
length(GEOgsea$pathway[GEOgsea$padj <= 0.05 & GEOgsea$ES > 0])
length(GEOgsea$pathway[GEOgsea$padj <= 0.05 & GEOgsea$ES < 0])
length(PPMIgsea$pathway[PPMIgsea$padj <= 0.05 & PPMIgsea$ES > 0])
length(PPMIgsea$pathway[PPMIgsea$padj <= 0.05 & PPMIgsea$ES < 0])

GEOmeta <- plotEnrichmentData(mito43[["REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES"]], GEOscore)
PPMImeta <- plotEnrichmentData(mito43[["REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES"]], PPMIscore)

with(GEOmeta,
    ggplot(data=curve) +
    geom_line(aes(x=rank, y=ES), color="green") +
    geom_ribbon(data=stats,
    mapping=aes(x=rank, ymin=0,
    ymax=stat/maxAbsStat*(spreadES/4)),
    fill="grey") +
    geom_segment(data=ticks,
    mapping=aes(x=rank, y=-spreadES/16,
    xend=rank, yend=spreadES/16),
    size=0.2) +
    geom_hline(yintercept=posES, colour="red", linetype="dashed") +
    geom_hline(yintercept=negES, colour="red", linetype="dashed") +
    geom_hline(yintercept=0, colour="black") +
    theme(
    panel.background = element_blank(),
    panel.grid.major=element_line(color="grey92")
    ) +
    labs(x="rank", y="enrichment score")))



########################



#Permutation testing
set.seed(01929)
totalgenes <- nrow(cisgenesgtex)
sample_size <- nrow(eitheror)
observed_true_count <- 9
num_permutations <- 10000
perm_true_counts <- numeric(num_permutations)
true_vector <- rep(TRUE, 140)
false_vector <- rep(FALSE, totalgenes - observed_true_count)
gene_vector <- c(true_vector, false_vector)


# Permutation test
for (i in 1:num_permutations) {
  # Shuffle the gene_vector to create a new permutation
  permuted_vector <- sample(gene_vector)
  
  # Count the number of TRUE values in the sample
  perm_true_count <- sum(permuted_vector[1:sample_size])
  perm_true_counts[i] <- perm_true_count
}

# Calculate p-value
p_value <- mean(perm_true_counts <= observed_true_count)

# Print the result
cat("Observed TRUE count:", observed_true_count, "\n")
cat("P-value:", p_value, "\n")

############################
#Summarise the eQTLs by the LD blocks

eqtls <- read.table("combinedeqtls.txt", sep = "\t", header = TRUE)
eqtlsummary <- eqtls %>%
  group_by(LDblock) %>%
  summarise(n = n(),
            n_mito = sum(impi))
write.csv(eqtlsummary, "combinedeqtlsummary.csv", row.names = FALSE, quote = FALSE)



#Creating venn diagram of all genes


sigGEO <- GEO[GEO$p_values2 <= 0.05,]
sigPPMI <- PPMI[PPMI$P.Value <= 0.05,]
sigGEO <- sigGEO[,1]
sigPPMI <- sigPPMI[,1]
combinedsig <- c(sigGEO, sigPPMI)


megalist <- Reduce(union, list(GEO$X, PPMI$X, counts$X, GTExexprnames$hgnc_symbol))
megalist <- as.data.frame(megalist)
write.table(megalist, "All possible background genes.txt", sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
megalist <- read.table("All possible background genes.txt", sep = "\t", header = FALSE)
colnames(megalist) <- c("genename")
megalist$cisgene <- megalist$genename %in% cisgenes$geneName
megalist$GTEx <- megalist$genename %in% GTExexprnames$hgnc_symbol
megalist$BRAINEAC <- megalist$genename %in% braineacexp$Name
megalist$GEO <- megalist$genename %in% GEO$X
megalist$PPMI <- megalist$genename %in% PPMI$X
megalist$eQTL <- megalist$genename %in% combined$gene
megalist$mito <- megalist$genename %in% impionly$Symbol
megalist$DEG <- megalist$genename %in% combinedsig
megalist$GTExeQTL <- megalist$genename %in% GTEx$gene
megalist$BRAINEACeQTL <- megalist$genename %in% BRAINEAC$gene
megalist$GEOsig <- megalist$genename %in% sigGEO
megalist$PPMIsig <- megalist$genename %in% sigPPMI
write.table(megalist, "All possible background genes annotated.txt", sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

megalist <- read.table("All possible background genes annotated.txt", sep = "\t", header = TRUE)

GTExnames <- megalist$genename[megalist$GTEx == TRUE]
BRAINEACnames <- megalist$genename[megalist$BRAINEAC == TRUE]
GEOnames <- megalist$genename[megalist$GEO == TRUE]
PPMInames <- megalist$genename[megalist$PPMI == TRUE]
cisgenenames <- megalist$genename[megalist$cisgene == TRUE]
mitonames <- megalist$genename[megalist$mito == TRUE]
DEGnames <- megalist$genename[megalist$DEG == TRUE]
eQTLnames <- megalist$genename[megalist$eQTL == TRUE]
GEOsignames <- megalist$genename[megalist$GEOsig == TRUE]
PPMIsignames <- megalist$genename[megalist$PPMIsig == TRUE]
GTExeQTLnames <- megalist$genename[megalist$GTExeQTL == TRUE]
BRAINEACeQTLnames <- megalist$genename[megalist$BRAINEACeQTL == TRUE]
GEOcissignames <- megalist$genename[megalist$GEOsig == TRUE & megalist$cisgene == TRUE]
PPMIcissignames <- megalist$genename[megalist$PPMIsig == TRUE & megalist$cisgene == TRUE]


fisherdata4 <- data.frame(
  "eQTL -" = c(length(megalist$genename[megalist$cisgene == TRUE & megalist$mito == FALSE & megalist$eQTL == FALSE]), length(megalist$genename[megalist$cisgene == TRUE & megalist$mito == TRUE & megalist$eQTL == FALSE])),
  "eQTL +" = c(length(megalist$genename[megalist$eQTL == TRUE & megalist$mito == FALSE]), length(megalist$genename[megalist$eQTL == TRUE & megalist$mito == TRUE])),
  row.names = c("NEMG -", "NEMG +"),
  stringsAsFactors = FALSE
)

fisher4 <- fisher.test(fisherdata4)
fisher4
########


eqtldegs <- megalist$genename[megalist$eQTL == TRUE & megalist$DEG == TRUE & megalist$mito == TRUE]
eqtldegs <- data.frame(eqtldegs)
eqtldegs$GEOZ <- GEO$Z[match(eqtldegs$eqtldegs, GEO$X)]
eqtldegs$PPMIZ <- PPMI$t[match(eqtldegs$eqtldegs, PPMI$X)]
eqtldegs$GEOP <- GEO$p_values2[match(eqtldegs$eqtldegs, GEO$X)]
eqtldegs$PPMIP <- PPMI$P.Value[match(eqtldegs$eqtldegs, PPMI$X)]
eqtldegs$LDblock <- combined$LDblock[match(eqtldegs$eqtldegs, combined$gene)]
combined <- combined[combined$impi == TRUE,]
combined <- combined[combined$gene %in% eqtldegs$eqtldegs,]
write.csv(eqtldegs, "eqtldegs.csv", row.names = FALSE, quote = FALSE)

######

upsetingput <- c(
  GTEx = (length(megalist$genename[megalist$GTExeQTL == TRUE]) + 1),
  BRAINEAC = length(megalist$genename[megalist$BRAINEACeQTL == TRUE]),
  GEO = length(megalist$genename[megalist$GEOsig == TRUE & megalist$cisgene == TRUE]),
  PPMI = length(megalist$genename[megalist$PPMIsig == TRUE & megalist$cisgene == TRUE]),
  eQTL = (length(megalist$genename[megalist$eQTL == TRUE]) + 1),
  DEG = length(megalist$genename[megalist$DEG == TRUE & megalist$cisgene == TRUE]),
  "GTEx&BRAINEAC" = length(megalist$genename[megalist$GTExeQTL == TRUE & megalist$BRAINEACeQTL == TRUE]),
  "GEO&PPMI" = length(megalist$genename[megalist$GEOsig == TRUE & megalist$PPMIsig == TRUE & megalist$cisgene == TRUE]),
  "eQTL&DEG" = length(megalist$genename[megalist$eQTL == TRUE & megalist$DEG == TRUE]),
  "GTEx&GEO" = length(megalist$genename[megalist$GTExeQTL == TRUE & megalist$GEOsig == TRUE]),
  "GTEx&PPMI" = length(megalist$genename[megalist$GTExeQTL == TRUE & megalist$PPMIsig == TRUE]),
  "BRAINEAC&DEG" = length(megalist$genename[megalist$BRAINEACeQTL == TRUE & megalist$DEG == TRUE]),
  "GTEx&DEG" = length(megalist$genename[megalist$GTExeQTL == TRUE & megalist$DEG == TRUE])
)
png("upsetplot.png", height = 800, width = 800)
upset(fromExpression(upsetingput), nsets = 6, sets = c("DEG", "eQTL", "GEO", "PPMI", "GTEx", "BRAINEAC"),
  keep.order = TRUE, order.by = "freq", decreasing = TRUE,
  sets.x.label = "Total interaction", mainbar.y.label = "Number of genes in set",
  set_size.show = FALSE, point.size = 6, line.size = 1.5, text.scale = 2,
  set_size.scale_max = 0
)
dev.off()
upset_data <- 
rownames(megalist) <- megalist$genename
megalist <- megalist[,-1]
lesserlist <- megalist[megalist$cisgene == TRUE,6:12]
lesserlist <- lesserlist[rowSums(lesserlist) >= 1, ]
categories <- colnames(lesserlist)
set_size(8, 3)
upset(lesserlist, categories)











#####
venndata <- list(
  GTEx = GTExnames,
  BRAINEAC = BRAINEACnames,
  GEO = GEOnames,
  PPMI = PPMInames
)

venndatasig <- list(
  GEO = GEOcissignames,
  PPMI = PPMIcissignames,
  GTEx = GTExeQTLnames,
  BRAINEAC = BRAINEACeQTLnames
)
png("vennallgenes.png")
ggplot() +
  geom_venn(venndata, textsize = 8) +
  theme_void() 
dev.off()



png("vennsigcisgenes.png")
ggplot() +
  geom_venn(venndatasig, textsize = 8) +
  theme_void()
dev.off()


install.packages("VennDiagram")
library(VennDiagram)

venn_plot <- venn.diagram(
  x = venndata,
  category.names = c(
,
  filename = "venn_diagram.png"
)

grid.draw(venn.plot)




########## block histograms
blockinformation <- read.csv("Block information.csv", header = TRUE)
plot <- ggplot(blockinformation, aes(x = widthkb, y = widthLD)) +
  geom_point(size = 0.5) +
  theme_bw() +
  xlab("Physical size (kb)") +
  ylab("Genetic length (LDU)")
png("blockscatter.png", width = 500, height = 500)
ggMarginal(plot, type = "histogram")
dev.off()

ggplot(blockinformation, aes(x=widthkb, y=widthLD)) +
  geom_point() +
  theme(legend.position = "none")


png("blockhistogram.png")
hist(blockinformation$widthkb, breaks = 10, main = "Histogram of block sizes", xlab = "Block size (kb)", ylab = "Frequency")
dev.off()
