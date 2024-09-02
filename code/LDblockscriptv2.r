#Setting working directory and loading packages
setwd("C:/Users/tb619/OneDrive - Imperial College London/PD_MSc_Project2024/Working materials/Visual Studio/Fixing blocks and making LDU0.1/")
install.packages("pak")
pak::pak("tidyverse/dplyr")
library(tidyverse)

#loading the table
ldustack<-read.csv("ldu_stacked.csv")
View(ldustack)
nall <- read.table("TableS2_v2.txt",head=T, sep="\t")

#assembling the hg37 data
lift_b35 <- as.data.frame(cbind(paste(rep("chr",dim(ldustack)[1]), ldu[,2], sep="")))
lift_b35[,2] <- as.integer(ldustack[,3]*1000)
lift_b35[,3] <- as.integer(ldustack[,3]*1000) + 1
lift_b35[,4] <- ldustack[,1]
colnames(lift_b35) <- c("chr","pos1", "pos2", "rsid")
head(lift_b35)
table(is.na(lift_b35[,2]))
write.table(lift_b35, file="lift_b35.bed", sep=" ", quote=F, row.names=F, col.names=F)
rm(lift_b35)

#the .bed file is then used in the ucsc liftover tool to convert to hg37 in a file named "b37lifted.bed". 52861 snps fail this process, but 2057626 are left.
b37<-read.table("b37lifted.bed")
b37<-subset(b37, select=-c(V3,V5))
colnames(b37)<-c("chromosome","pos37","snp")

#import the hg37 values into the ldu file
ldustack<-merge(ldustack,b37, by = "snp")
rm(b37)

#insert the missing snps from nall into ldu data frames
nall$inldu<-nall[,1] %in% ldustack[,1]
importing<-subset(nall, nall$inldu==FALSE)
importing<-subset(importing, select=c(snp,chr,BP_b37))
colnames(importing)<-c("snp","chr","pos37")
class(importing$pos37)
importing$pos37<-gsub(",","",importing$pos37)
importing$pos37<-as.integer(importing$pos37)
class(importing$pos37)
ldustack<-bind_rows(ldustack,importing)
ldustack$pos37<-ldustack$pos37/1000
rm(importing)

#fix the ordering
fixorder <- order(ldustack$chr, ldustack$pos37)
ldustack<-ldustack[fixorder,]
rownames(ldustack)<- NULL
ldustack<-subset(ldustack, select=-chromosome)
rm(fixorder)

#mark the lead snps
ldustack$lead<- ldustack[,1] %in% nall[,1]
table(ldustack$lead)
ldustack[ldustack$lead==T,]

#give the snps missing LDU a value
for (i in 1:(nrow(ldustack) - 1)) {
  if (is.na(ldustack[i, "ldu"])) {
    ldustack[i, "ldu"] <- mean(ldustack[c(i-1, i+1), "ldu"], na.rm = TRUE)
  }
}
rm(i)
#select all SNPs which have the exact same LDU as the lead SNPs
leadrows <- subset(ldustack, lead == TRUE)
result <- list()
for (i in 1:nrow(leadrows)) {
  currentchr <- leadrows[i, "chr"]
  matchingrows <- subset(ldustack, chr == currentchr & leadrows[i, "ldu"] == ldu)
  result[[i]] <- matchingrows
}
result<-unique(result)
ldu0 <- do.call(rbind, result)

#check there is the right number of lead SNPs still
sum(ldu0$lead==TRUE)
#export the file and clean up
write.csv(ldu0, file="LDU0 around lead SNPs.csv", row.names = FALSE)
rm(lead_rows,matching_rows,result,current_chr,i)

#next, a file with all SNPs with ±0.5LDU
library(dplyr)

leadrows <- filter(ldustack, lead == TRUE)
ldu0.5 <- ldustack %>%
    filter(chr %in% leadrows$chr, abs(ldu - leadrows$ldu) <= 0.5) %>%
    unique()

lead_count <- sum(ldu0.5$lead == TRUE)
write.csv(ldu0.5,file="LDU0.5 around lead SNPs.csv", row.names = FALSE)
rm(leadrows,matchingrows,result,currentchr,i)

#group SNps by LD block
ldu0<-read.csv("LDU0 around lead SNPs.csv")
ldu0.5<-read.csv("LDU0.5 around lead SNPs.csv")
ldu0 <- mutate(ldu0, LDblock = group_indices(ldu0, chr, ldu))


#assigning the SNPs to a block was done manually, as there is messy boundaries between blocks
ldu0.5$LDblock<-NA
ldu0.5[1:72,"LDblock"]<-1
ldu0.5[73:110,"LDblock"]<-2
ldu0.5[111:206,"LDblock"]<-3
ldu0.5[207:260,"LDblock"]<-4
ldu0.5[261:343,"LDblock"]<-5
ldu0.5[344:421,"LDblock"]<-6
ldu0.5[422:456,"LDblock"]<-7
ldu0.5[457:485,"LDblock"]<-8
ldu0.5[486:505,"LDblock"]<-9
ldu0.5[506:662,"LDblock"]<-10
ldu0.5[663:797,"LDblock"]<-11
ldu0.5[798:883,"LDblock"]<-12
ldu0.5[884:919,"LDblock"]<-13
ldu0.5[920:995,"LDblock"]<-14
ldu0.5[996:1048,"LDblock"]<-15
ldu0.5[1049:1331,"LDblock"]<-16
ldu0.5[1332:1447,"LDblock"]<-17
ldu0.5[1448:1528,"LDblock"]<-18
ldu0.5[1529:1814,"LDblock"]<-19
ldu0.5[1815:1913,"LDblock"]<-20
#Two leads in this block
ldu0.5[1914:1994,"LDblock"]<-21
ldu0.5[1995:2033,"LDblock"]<-22
ldu0.5[2034:2166,"LDblock"]<-23
ldu0.5[2167:2278,"LDblock"]<-24
ldu0.5[2279:2304,"LDblock"]<-25
#these two need to be careful
ldu0.5[2305:2411,"LDblock"]<-26
#Two leads in this block
ldu0.5[2412:2600,"LDblock"]<-27
ldu0.5[2601:2680,"LDblock"]<-28
ldu0.5[2681:2743,"LDblock"]<-29
ldu0.5[2744:3034,"LDblock"]<-30
ldu0.5[3035:3229,"LDblock"]<-31
ldu0.5[3230:3299,"LDblock"]<-32
ldu0.5[3300:3654,"LDblock"]<-33
ldu0.5[3655:4066,"LDblock"]<-34
ldu0.5[4067:4261,"LDblock"]<-35
ldu0.5[4262:4382,"LDblock"]<-36
ldu0.5[4383:4465,"LDblock"]<-37
ldu0.5[4466:4651,"LDblock"]<-38
ldu0.5[4652:4801,"LDblock"]<-39
ldu0.5[4802:5193,"LDblock"]<-40
ldu0.5[5194:5249,"LDblock"]<-41
ldu0.5[5250:5332,"LDblock"]<-42
ldu0.5[5333:5395,"LDblock"]<-43
ldu0.5[5396:5458,"LDblock"]<-44
ldu0.5[5459:5537,"LDblock"]<-45
ldu0.5[5538:5657,"LDblock"]<-46
ldu0.5[5658:5832,"LDblock"]<-47
ldu0.5[5833:5936,"LDblock"]<-48
ldu0.5[5937:6021,"LDblock"]<-49
#these two need care
ldu0.5[6022:6074,"LDblock"]<-50
ldu0.5[6075:6113,"LDblock"]<-51
ldu0.5[6114:6220,"LDblock"]<-52
ldu0.5[6221:6299,"LDblock"]<-53
ldu0.5[6300:6319,"LDblock"]<-54
ldu0.5[6320:6476,"LDblock"]<-55
ldu0.5[6477:6795,"LDblock"]<-56
ldu0.5[6796:6910,"LDblock"]<-57
ldu0.5[6911:6975,"LDblock"]<-58
ldu0.5[6976:7019,"LDblock"]<-59
ldu0.5[7020:7040,"LDblock"]<-60
ldu0.5[7041:7112,"LDblock"]<-61
ldu0.5[7113:7342,"LDblock"]<-62
ldu0.5[7343:7434,"LDblock"]<-63
ldu0.5[7435:7584,"LDblock"]<-64
ldu0.5[7585:7788,"LDblock"]<-65
ldu0.5[7789:7830,"LDblock"]<-66
ldu0.5[7831:7859,"LDblock"]<-67
ldu0.5[7860:7935,"LDblock"]<-68
ldu0.5[7936:8019,"LDblock"]<-69
ldu0.5[8020:8101,"LDblock"]<-70
ldu0.5[8102:8167,"LDblock"]<-71
ldu0.5[8168:8211,"LDblock"]<-72
ldu0.5[8212:8267,"LDblock"]<-73
ldu0.5[8268:8328,"LDblock"]<-74
ldu0.5[8329:8364,"LDblock"]<-75
ldu0.5[8365:8371,"LDblock"]<-76
ldu0.5[8372:9235,"LDblock"]<-77
ldu0.5[9236:9261,"LDblock"]<-78
ldu0.5[9262:9346,"LDblock"]<-79
ldu0.5[9347:9372,"LDblock"]<-80
ldu0.5[9373:9571,"LDblock"]<-81
ldu0.5[9572:9729,"LDblock"]<-82
ldu0.5[9730:9805,"LDblock"]<-83
ldu0.5[9806:9822,"LDblock"]<-84
ldu0.5[9823:9862,"LDblock"]<-85
ldu0.5[9863:10027,"LDblock"]<-86
write.csv(ldu0, file="LDU0 around lead SNPs.csv", row.names = FALSE)
write.csv(ldu0.5, file="LDU0.5 around lead SNPs.csv", row.names = FALSE)

#Create summary statistics for the two data frames

sumtableldu0 <- ldu0 %>%
  group_by(LDblock) %>%
  summarize(
    number = n(),
    maxLD = max(ldu),
    minLD = min(ldu),
    medLD = median(ldu),
    widthLD = maxLD - minLD,
    maxkb = max(pos37),
    minkb = min(pos37),
    medkb = median(pos37),
    widthkb = maxkb - minkb
  )

write.csv(sumtableldu0, file = "Summary for LDU0.csv", row.names = FALSE)

sumtableldu0.5 <- ldu0.5 %>%
  group_by(LDblock) %>%
  summarize(
    number = n(),
    maxLD = max(ldu),
    minLD = min(ldu),
    medLD = median(ldu),
    widthLD = maxLD - minLD,
    maxkb = max(pos37),
    minkb = min(pos37),
    medkb = median(pos37),
    widthkb = maxkb - minkb
  )
write.csv(sumtableldu0.5, file="Summary for LDU0.5.csv", row.names = FALSE)


plot(sumtableldu0.5$widthLD, sumtableldu0.5$widthkb, cex=0.8, pch=16, main = "LD block width in LDU vs kb", ylab = "Width in kb", xlab="Width in LDU")
abline(lm(sumtableldu0.5$widthkb ~ sumtableldu0.5$widthLD))
plot(sumtableldu0.5$number, sumtableldu0.5$widthkb, cex=0.8, pch=16, main="LD block width vs group size", xlab = "Number of SNPs in block", ylab = "Block width in kb")
abline(lm(sumtableldu0.5$widthkb~sumtableldu0.5$number))
library(rgl)
plot3d(sumtableldu0.5$widthkb, sumtableldu0.5$widthLD, sumtableldu0.5$number, size=1, xlab="Width in kb", ylab="Width in LDU", zlab="Number of SNPs in block", type='s')
htmlwidgets::saveWidget(rglwidget(width = 800, height = 800),
                        file="3D scatter for LDU0.5.html",
                        libdir="libs",
                        selfcontained = TRUE)
boxplot(sumtableldu0$widthkb)
boxplot(sumtableldu0.5$widthkb)
