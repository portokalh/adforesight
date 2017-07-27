
setwd('/Users/alex/AlexBadea_mycodes/rplay')
manganese_t1 <- read.csv(file="All_BehaviorChiMn2.csv",head=TRUE,sep=",")
manganese_t1 <- read.csv(file="All_BehaviorChiMn.csv",head=TRUE,sep=",")

#manganese_t1_cvn<-manganese_t1$MnOriMean[manganese_t1$$genotype==1]
#manganese_t1_nos<-manganese_t1$MnOriMean[manganese_t1$$genotype==0]

#T1w<-manganese_t1$MnOriMean
#genotype<-manganese_t1$$genotype

T1wMnGenotype_t<-t.test(MnOriMean~genotype,data=manganese_t1,alternative="two.sided")
QSMGlobal_t<-t.test(chi~genotype,data=manganese_t1)
T1wMnSex_t<-t.test(MnOriMean~sex,data=manganese_t1)
QSMSex_t<-t.test(chi~sex,data=manganese_t1)

pdf(file='T1wSignalGlobalByGenotype', width=4, height=4, pointsize=12)
boxplot(MnOriMean~genotype,data=manganese_t1, main="T1w Signal Brain", col=(c("white","red")),
        xlab="Genotype", ylab="T1w Signal Intensity")
dev.off()


pdf(file='QSMGlobalByGenotype', width=4, height=4, pointsize=12)
boxplot(chi~genotype,data=manganese_t1, main="QSM Brain", col=(c("white","red")),
        xlab="Genotype", ylab="QSM")
dev.off()

pdf(file='T1wSignalGlobalBySex', width=4, height=4, pointsize=12)
boxplot(MnOriMean~sex,data=manganese_t1, main="T1w Signal Brain", col=(c("white","red")),
        xlab="Sex", ylab="T1w Signal Intensity")
dev.off()

pdf(file='QSMGlobalBySex', width=4, height=4, pointsize=12)
boxplot(chi~sex,data=manganese_t1, main="QSM Brain", col=(c("white","red")),
        xlab="Sex", ylab="QSM")
dev.off()


