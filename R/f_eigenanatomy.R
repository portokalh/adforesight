## used for Fig4; here for T1w intensity contrast
## load the other contrasts for the different inputs (def, qsm)
library( knitr )
library(ANTsR)

setwd( '/Users/omega/alex/adforesight' )
output.path <- '/Users/omega/alex/adforesight/mydata/outdata/sd2_projall_noscale/'
if (dir.exists(output.path)){ 1} else {dir.create(output.path, recursive=TRUE)}

# Load in Behavior and Imaging Data
behavior <- read.csv('./mydata/All_Behavior.csv')
labled.set <-read.csv('./mydata/legendsCHASS2symmetric.csv')
labeled.brain.img <-  antsImageRead('./mydata/MDT_labels_chass_symmetric.nii.gz')
mask <-antsImageRead('./mydata/MDT_mask_e3.nii')
mang_files<-list.files(path = "./mydata/imdata/", pattern = "T2_to_MDT",full.names = T,recursive = T)
mang_mat <-imagesToMatrix(mang_files,mask)

rows.test <-as.integer(c(2,3,4,9,15,17))
rows.train <-as.integer(c(23,8,18,13,14,19,16,12,11,24,20,7,6,1,22,10,5,21))

mang.train <- mang_mat[rows.train, ]
mang.test <- mang_mat[rows.test, ]
behav.train <- behavior[rows.train, ]
behav.test <-behavior[rows.test, ]

eanat_mang_un<-sparseDecom( inmatrix=mang.train, inmask=mask, nvecs=50,
                            sparseness=0.005, cthresh=250, its=5, mycoption=0 )
jeanat_mang<-joinEigenanatomy(mang.train , mask, eanat_mang_un$eigenanatomyimages,
                              0.1, joinMethod='multilevel' )
useeig_mang<-jeanat_mang$fusedlist
avgmat_mang<-abs(imageListToMatrix( useeig_mang , mask ))
avgmat_mang<-avgmat_mang/rowSums(abs(avgmat_mang))
imgmat_train_mang<-(mang.train %*% t(avgmat_mang)  )
imgmat_test_mang<-(mang.test %*% t(avgmat_mang)  )

dist4.train <-behav.train[,'d4']
dist4.test <-behav.test[,'d4']
projs.train <- data.frame(cbind(dist4.train,imgmat_train_mang)) # column combine the behavior wth the projections
colnames(projs.train) <- c('Dist_4', paste0('Proj', c(1:ncol(imgmat_train_mang)))) # insert column names
projs.test <- data.frame(cbind(dist4.test,imgmat_test_mang)) # column combine the behavior wth the projections
colnames(projs.test) <- c('Dist_4', paste0('Proj', c(1:ncol(imgmat_test_mang)))) # insert column names
mylm <- lm('Dist_4 ~ .', data=projs.train) # behavior correlation with the number of projections
distpred <- predict.lm(mylm, newdata=projs.test) # based on the linear model predict the distances for the same day
modsum <-summary(lm(distpred~dist4.test))
r2 <- modsum$adj.r.squared
my.p <- modsum$coefficients[2,4]
plot(dist4.test, distpred, xlab = 'Real Swim Distance on Day 4', ylab = 'Predicted Swim Distance on day 4', 
     main='Predicted vs. Real Swim Distance on Day 4', ylim = c(0,1000),xlim = c(0,1000)) # generate plot
abline(lm(distpred~dist4.test))
abline(a = 0, b =1, pch=22, lty=2, col="red")
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE),
                   list(MYVALUE = format(r2, digits = 3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE),
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
legend('topleft', legend = rp, bty = 'n')



mylm <- lm('Dist_4 ~ Proj6', data=projs.train) # behavior correlation with the number of projections
distpred <- predict.lm(mylm, newdata=projs.test) # based on the linear model predict the distances for the same day
modsum <-summary(lm(distpred~dist4.test))
r2 <- modsum$adj.r.squared
my.p <- modsum$coefficients[2,4]
plot(dist4.test, distpred, xlab = 'Real Swim Distance on Day 4', ylab = 'Predicted Swim Distance on day 4', 
     main='Predicted vs. Real Swim Distance on Day 4 Eigregion 6', ylim = c(0,1000),xlim = c(0,1000)) # generate plot
abline(lm(distpred~dist4.test))
abline(a = 0, b =1, pch=22, lty=2, col="red")
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE),
                   list(MYVALUE = format(r2, digits = 3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE),
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
legend('topleft', legend = rp, bty = 'n')