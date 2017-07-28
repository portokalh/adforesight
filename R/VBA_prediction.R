library( knitr )
library(ANTsR)

# Load in Behavior and Imaging Data
behavior <- read.csv("/Users/omega/Documents/Natalie/SCCAN_tutorial/Eig_zscore/data/All_Behavior.csv")
mask <-antsImageRead('/Users/omega/Documents/Natalie/SCCAN_tutorial/Eig_zscore/data/MDT_mask_e3.nii')
jac.result <- antsImageRead('/Users/omega/Documents/Natalie/VBA/spm_results/JAC/jac_result_sub.nii.gz')
pos.clust.jac <-thresholdImage(jac.result,2,17)
neg.clust.jac <-thresholdImage(jac.result,-15,-2)
#plot.antsImage(pos.clust.jac)
#plot.antsImage(neg.clust.jac)
#antsImageWrite(pos.clust.jac, '/Users/omega/Documents/Natalie/VBA/spm_results/JAC/jac_pos_result.nii.gz')
#antsImageWrite(neg.clust.jac, '/Users/omega/Documents/Natalie/VBA/spm_results/JAC/jac_neg_result.nii.gz')
jac_files<-list.files(path = "/Users/omega/Documents/Natalie/SCCAN_tutorial/Eig_zscore/data", pattern = "jac_to_MDT",full.names = T,recursive = T)
jac_mat <-imagesToMatrix(jac_files,mask)
back <-matrixToImages(jac_mat,mask)

rows.test <-as.integer(c(2,3,4,9,15,17))
rows.train <-as.integer(c(23,8,18,13,14,19,16,12,11,24,20,7,6,1,22,10,5,21))

jac.train <- jac_mat[rows.train, ]
jac.test <- jac_mat [rows.test, ]
behav.train <- behavior[rows.train, ]
behav.test <-behavior[rows.test, ]
initmat <- matrix(0,nrow = 1, ncol = 492797)

#b <- back[[1]]
vec<-( jac.result[ mask == 1 ] >= 2 )
initmat[1,]<-as.numeric( vec )

kp <-matrix(0, nrow = length(rows.train), ncol = 492797)
for (i in 1:length(rows.train)){
kp[i,] <- initmat * jac.train[i, ]
}

imgmat_jac<-as.matrix(colMeans(kp))
imgpredtrain_jac<-jac.train %*% (imgmat_jac)
imgpredtest_jac<-jac.test %*% (imgmat_jac)

dist4.train <-behav.train[,'d4']
dist4.test <-behav.test[,'d4']
projs.train <- data.frame(cbind(dist4.train,imgpredtrain_jac)) # column combine the behavior wth the projections
colnames(projs.train) <- c('Dist_4', paste0('Proj', c(1:ncol(imgpredtrain_jac)))) # insert column names
projs.test <- data.frame(cbind(dist4.test,imgpredtest_jac)) # column combine the behavior wth the projections
colnames(projs.test) <- c('Dist_4', paste0('Proj', c(1:ncol(imgpredtest_jac)))) # insert column names
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

###### SUSCEPT #####
suscept.result <- antsImageRead('/Users/omega/Documents/Natalie/VBA/spm_results/Chi_subtracted/suscept_result_sub.nii.gz')
# pos.clust.suscept <-thresholdImage(suscept.result,2,10)
# neg.clust.suscept <-thresholdImage(suscept.result,-10,-2)
# plot.antsImage(pos.clust.suscept)
# plot.antsImage(neg.clust.suscept)
# antsImageWrite(pos.clust.suscept, '/Users/omega/Documents/Natalie/VBA/spm_results/Chi_subtracted/suscept_pos_result.nii.gz')
# antsImageWrite(neg.clust.suscept, '/Users/omega/Documents/Natalie/VBA/spm_results/Chi_subtracted/suscept_neg_result.nii.gz')
suscept_files<-list.files(path = "/Users/omega/Documents/Natalie/SCCAN_tutorial/Eig_zscore/data", pattern = "X_to_MDT",full.names = T,recursive = T)
suscept_mat <-imagesToMatrix(suscept_files,mask)

rows.test <-as.integer(c(2,3,4,9,15,17))
rows.train <-as.integer(c(23,8,18,13,14,19,16,12,11,24,20,7,6,1,22,10,5,21))

suscept.train <- suscept_mat[rows.train, ]
suscept.test <- suscept_mat [rows.test, ]
behav.train <- behavior[rows.train, ]
behav.test <-behavior[rows.test, ]
initmat <- matrix(0,nrow = 1, ncol = 492797)

vec<-( suscept.result[ mask == 1 ] <= -2 )
initmat[1,]<-as.numeric( vec )

kp <-matrix(0, nrow = length(rows.train), ncol = 492797)
for (i in 1:length(rows.train)){
  kp[i,] <- initmat * suscept.train[i, ]
}

imgmat_suscept<-as.matrix(colMeans(kp))
imgpredtrain_suscept<-suscept.train %*% (imgmat_suscept)
imgpredtest_suscept<-suscept.test %*% (imgmat_suscept)

dist4.train <-behav.train[,'d4']
dist4.test <-behav.test[,'d4']
projs.train <- data.frame(cbind(dist4.train,imgpredtrain_suscept)) # column combine the behavior wth the projections
colnames(projs.train) <- c('Dist_4', paste0('Proj', c(1:ncol(imgpredtrain_suscept)))) # insert column names
projs.test <- data.frame(cbind(dist4.test,imgpredtest_suscept)) # column combine the behavior wth the projections
colnames(projs.test) <- c('Dist_4', paste0('Proj', c(1:ncol(imgpredtest_suscept)))) # insert column names
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

#### T1-Mang ####
mang.result <- antsImageRead('/Users/omega/Documents/Natalie/VBA/spm_results/T1/mang_result_sub.nii.gz')
pos.clust.mang <-thresholdImage(mang.result,2,10)
neg.clust.mang <-thresholdImage(mang.result,-10,-2)
plot.antsImage(pos.clust.mang)
plot.antsImage(neg.clust.mang)
antsImageWrite(pos.clust.mang, '/Users/omega/Documents/Natalie/VBA/spm_results/T1/mang_pos_result.nii.gz')
antsImageWrite(neg.clust.mang, '/Users/omega/Documents/Natalie/VBA/spm_results/T1/mang_neg_result.nii.gz')
mang_files<-list.files(path = "/Users/omega/Documents/Natalie/SCCAN_tutorial/Eig_zscore/data", pattern = "T2_to_MDT",full.names = T,recursive = T)
mang_mat <-imagesToMatrix(mang_files,mask)

rows.test <-as.integer(c(2,3,4,9,15,17))
rows.train <-as.integer(c(23,8,18,13,14,19,16,12,11,24,20,7,6,1,22,10,5,21))

mang.train <- mang_mat[rows.train, ]
mang.test <- mang_mat [rows.test, ]
behav.train <- behavior[rows.train, ]
behav.test <-behavior[rows.test, ]
initmat <- matrix(0,nrow = 1, ncol = 492797)

vec<-( mang.result[ mask == 1 ] >= 2 )
initmat[1,]<-as.numeric( vec )

kp <-matrix(0, nrow = length(rows.train), ncol = 492797)
for (i in 1:length(rows.train)){
  kp[i,] <- initmat * mang.train[i, ]
}

imgmat_mang<-as.matrix(colMeans(kp))
imgpredtrain_mang<-mang.train %*% (imgmat_mang)
imgpredtest_mang<-mang.test %*% (imgmat_mang)

dist4.train <-behav.train[,'d4']
dist4.test <-behav.test[,'d4']
projs.train <- data.frame(cbind(dist4.train,imgpredtrain_mang)) # column combine the behavior wth the projections
colnames(projs.train) <- c('Dist_4', paste0('Proj', c(1:ncol(imgpredtrain_mang)))) # insert column names
projs.test <- data.frame(cbind(dist4.test,imgpredtest_mang)) # column combine the behavior wth the projections
colnames(projs.test) <- c('Dist_4', paste0('Proj', c(1:ncol(imgpredtest_mang)))) # insert column names
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