library(broom)
library(gridExtra)
library(pander)


#setwd('/Users/alex/AlexBadea_mycodes/rplay')

# set working directory here
setwd('~/GitHub/adforesight')


#this one has logical vars for sex and genotype
manganese_t1 <- read.csv(file="./R/All_BehaviorChiMn2.csv",head=TRUE,sep=",")


#this old one has numeric not logical vars for sex and genotype
#manganese_t1 <- read.csv(file="All_BehaviorChiMn.csv",head=TRUE,sep=",")

#manganese_t1_cvn<-manganese_t1$MnOriMean[manganese_t1$$genotype==1]
#manganese_t1_nos<-manganese_t1$MnOriMean[manganese_t1$$genotype==0]

#T1w<-manganese_t1$MnOriMean
#genotype<-manganese_t1$$genotype
#remove first arg for cbind
VolGlobalGenotype_t<-cbind("VolGenotype",tidy(t.test(BrainVol~genotype,data=manganese_t1)))
VolGlobalSex_t<-cbind("VolSex",tidy(t.test(BrainVol~sex,data=manganese_t1)))
T1wMnGenotype_t<-cbind("T1wGenotype",tidy(t.test(MnOriMean~genotype,data=manganese_t1,alternative="two.sided")))
T1wMnNormedGlobalGenotype_t<-cbind("T1wNSex",tidy(t.test(MnNormedMean~genotype,data=manganese_t1,alternative="two.sided")))
QSMGenotypeGlobal_t<-cbind("QSMGenotype",tidy(t.test(BrainChi~genotype,data=manganese_t1)))
T1wMnGlobalSex_t<-cbind("T1wSex",tidy(t.test(MnOriMean~sex,data=manganese_t1)))
T1wMnNormedGlobalSex_t<-cbind("T1wNSex",tidy(t.test(MnNormedMean~sex,data=manganese_t1)))
QSMGlobalSex_t<-cbind("QSMSex",tidy(t.test(BrainChi~sex,data=manganese_t1)))

#write global stats to file
allglobals <- rbind(VolGlobalGenotype_t,VolGlobalSex_t,T1wMnNormedGlobalGenotype_t,T1wMnNormedGlobalSex_t,QSMGenotypeGlobal_t,QSMGlobalSex_t)
rownames(allglobals) <- c("VolGlobalGenotype","VolGlobalSex", "T1wGlobalGenotype","T1wGlobalSex","QSMGlobalGenotype","QSMGlobalSex")
colnames(allglobals) <- c( "Obs",  "MeanDiff", "Mean1",	"Mean2",	"statistic",	"p.value", 	"df",	"conf.low",	"conf.high",	"method",	"alternative")

output.globalstatsfile <- file.path('./ROISTATS/global_stats.csv')
output.file<-write.table(allglobals, file = output.globalstatsfile , append = FALSE, sep = ",", col.names = TRUE, row.names = c("VolGenotype","VolSex", "T1wGenotype","T1wSex","QSMGenotype","QSMSex"))
#output.file<-write.table(allglobals, file = output.globalstatsfile , append = FALSE, sep = ",", col.names = c("Obs", "MeanDiff", "Mean1",	"Mean2",	"statistic",	"p.value", 	"df",	"conf.low",	"conf.high",	"method",	"alternative"), row.names = TRUE)
# concatenate my stats and write to csv; shifts my column titles to the left dag nabit

pandoc.table(allglobals[1:7, 1:11],caption = "Global Stats"))
knitr::kable(head(allglobals))


pdf(file='BrainVolGlobalByGenotype.pdf', width=4, height=4, pointsize=12)
boxplot(BrainVol~genotype,data=manganese_t1, main="Brain Vol (mm3)", col=(c("white","red")),
        xlab="Genotype", ylab="Brain Vol (mm3)")
dev.off()

pdf(file='BrainT1wGlobalByGenotype.pdf', width=4, height=4, pointsize=12)
boxplot(MnOriMean~genotype,data=manganese_t1, main="T1w Signal Brain (AU)", col=(c("white","red")),
        xlab="Genotype", ylab="T1w Signal Intensity")
dev.off()


pdf(file='BrainQSMGlobalByGenotype.pdf', width=4, height=4, pointsize=12)
boxplot(BrainChi~genotype,data=manganese_t1, main="QSM Brain (ppm)", col=(c("white","red")),
        xlab="Genotype", ylab="QSM")
dev.off()

pdf(file='BrainVolGlobalGlobalBySex.pdf', width=4, height=4, pointsize=12)
boxplot(BrainVol~sex,data=manganese_t1, main="Brain Vol (mm3)", col=(c("white","red")),
        xlab="Sex", ylab="QSM")
dev.off()

pdf(file='BrainT1wGlobalBySex.pdf', width=4, height=4, pointsize=12)
boxplot(MnOriMean~sex,data=manganese_t1, main="T1w Signal Brain (AU)", col=(c("white","red")),
        xlab="Sex", ylab="T1w Signal Intensity")
dev.off()

pdf(file='BrainQSMGlobalGlobalBySex.pdf', width=4, height=4, pointsize=12)
boxplot(BrainChi~sex,data=manganese_t1, main="QSM Brain (ppm)", col=(c("white","red")),
        xlab="Sex", ylab="QSM")
dev.off()

#locations for ROI-s
#'/Users/alex/brain_data/MEMRI_CVN_NOS/MEMRI_labels/labels_CVN_N0S2/'
#%interesting ROIs
#labelL.ind <- c(19,20,42,47,51,54,57,59,62,63,64,65,73,81,83,91,118,120,122)
#labelR.ind<-labelL.ind+1000