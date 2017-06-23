library(ANTsR)
library(knitr)

setwd( '/Users/alex/GitHub/adforesight' )
output.path <- '/Users/alex/GitHub/adforesight/mydata/outdata/sd2_projall_noscale/'
#setwd( '/Users/omega/alex/adforesight' )
#output.path <- '/Users/omega/alex/adforesight/mydata/outdata/sd2_projall_noscale/'


behavior <- read.csv("/Users/omega/Documents/Natalie/SCCAN_tutorial/Eig_zscore/data/All_Behavior.csv")
adind <-c(1,2,8,9,10,12,13,14,15,16,19)
conind<-c(3,4,5,6,7,11,17,18,20,21,22,23,24)
ad_behavior <-behavior[adind,]
con_behavior <-behavior[conind,]

#figure 1 Swim Distance, Swim time, swim velocity
a <- c(0.6,1.6,2.6,3.6,4.6)
b<-a-.1
png(filename = '/Users/omega/Documents/Natalie/Thesis/Figures/figures_boxplots/Swim_Distance_ConvsNOS2_slope.png')
ad_dist <-ad_behavior[, c(2:6)]
cd_dist <-con_behavior[, c(2:6)]
time <-1:5
one<-boxplot(ad_dist, ylab = 'Swim Distance (cm)', border="red", boxwex = 0.35 ,outline = T, ylim = c(0,1400))
med_addist1 <- median(ad_dist$d1)
med_addist2 <- median(ad_dist$d2)
med_addist3 <- median(ad_dist$d3)
med_addist4 <- median(ad_dist$d4)
med_addist5 <- median(ad_dist$d5)
med_addists<-c(med_addist1,med_addist2,med_addist3,med_addist4,med_addist5)
abline(lm(med_addists~time), lty =  1 , lwd = 2 , col = 'red')
two<-boxplot(cd_dist, ylab = 'Swim Distance (cm)', main = "Average Swim Distance Traveled", border="black", xaxt='n', boxwex = 0.35 , add = T, at = a, outline = T)
med_cddist1 <- median(cd_dist$d1)
med_cddist2 <- median(cd_dist$d2)
med_cddist3 <- median(cd_dist$d3)
med_cddist4 <- median(cd_dist$d4)
med_cddist5 <- median(cd_dist$d5)
med_cddists<-c(med_cddist1,med_cddist2,med_cddist3,med_cddist4,med_cddist5)
abline(lm(med_cddists~time), lty =  1 , lwd = 2 , col = 'black')
legend("bottomleft", title = 'Genotype', legend = c('NOS2', 'CVN'), fill=c('black','red'), horiz=TRUE)
dev.off()


png(filename = '/Users/omega/Documents/Natalie/Thesis/Figures/figures_boxplots/Swim_Time_ConvsNOS2_slope.png')
ad_time <-ad_behavior[, c(7:11)]
cd_time <-con_behavior[, c(7:11)]
boxplot(ad_time, ylab = 'Swim Time (sec)', border="red", boxwex = 0.35, outline = T, ylim= c(0,70))
med_adtime1 <- median(ad_time$s1)
med_adtime2 <- median(ad_time$s2)
med_adtime3 <- median(ad_time$s3)
med_adtime4 <- median(ad_time$s4)
med_adtime5 <- median(ad_time$s5)
med_adtimes<-c(med_adtime1,med_adtime2,med_adtime3,med_adtime4,med_adtime5)
abline(lm(med_adtimes~time), lty =  1 , lwd = 2 , col = 'red')
boxplot(cd_time, ylab = 'Swim Time (sec)', main = "Average Swim Time",  border="black",xaxt='n', boxwex = 0.35 , add = T, at = b, outline = T)
med_cdtime1 <- median(cd_time$s1)
med_cdtime2 <- median(cd_time$s2)
med_cdtime3 <- median(cd_time$s3)
med_cdtime4 <- median(cd_time$s4)
med_cdtime5 <- median(cd_time$s5)
med_cdtimes<-c(med_cdtime1,med_cdtime2,med_cdtime3,med_cdtime4,med_cdtime5)
abline(lm(med_cdtimes~time), lty =  1 , lwd = 2 , col = 'black')
legend("bottomleft", title = 'Genotype', legend = c('NOS2', 'CVN'), fill=c('black','red'), horiz=TRUE)
dev.off()

png(filename = '/Users/omega/Documents/Natalie/Thesis/Figures/figures_boxplots/Swim_Velocity_ConvsNOS2_slope.png')
ad_vel <-ad_behavior[, c(12:16)]
cd_vel <-con_behavior[, c(12:16)]
boxplot(ad_vel, ylab = 'Swim Velocity (cm/sec)', border="red", boxwex = 0.35, outline = T,ylim = c(10,25))
med_advel1 <- median(ad_vel$v1)
med_advel2 <- median(ad_vel$v2)
med_advel3 <- median(ad_vel$v3)
med_advel4 <- median(ad_vel$v4)
med_advel5 <- median(ad_vel$v5)
med_advels<-c(med_advel1,med_advel2,med_advel3,med_advel4,med_advel5)
abline(lm(med_advels~time), lty =  1 , lwd = 2 , col = 'red')
boxplot(cd_vel, ylab = 'Swim Velocity (cm/sec)', main = "Average Swim Velocity", border="black", xaxt='n', boxwex = 0.35 , add = T, at = a, outline = T)
med_cdvel1 <- median(cd_vel$v1)
med_cdvel2 <- median(cd_vel$v2)
med_cdvel3 <- median(cd_vel$v3)
med_cdvel4 <- median(cd_vel$v4)
med_cdvel5 <- median(cd_vel$v5)
med_cdvels<-c(med_cdvel1,med_cdvel2,med_cdvel3,med_cdvel4,med_cdvel5)
abline(lm(med_cdvels~time), lty =  1 , lwd = 2 , col = 'black')
legend("bottomleft", title = 'Genotype', legend = c('NOS2', 'CVN'), fill=c('black','red'), horiz=TRUE)
dev.off()