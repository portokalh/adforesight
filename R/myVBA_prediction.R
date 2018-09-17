library( knitr )
library(ANTsR)
library(ANTsRCore)
library(caret)

# Load in Behavior and Imaging Data
mypath<-'/Users/alex/GitHub/adforesight/' #'/Volumes/CivmUsers/omega/alex/GitHub/adforesight/'
# Load in Behavior and Imaging Data
#setwd( '/Users/alex/GitHub/adforesight/' )
setwd(mypath)
extension<-'VBA_Mn_pos_pred'
output.path <- paste(mypath,'/mydata/outdata/',extension, '/', sep='') #sd2_projall_noscale/'
vbainput.path<-'/Users/alex/brain_data/MEMRI_CVN_NOS/MEMRI_VBA/'
#setwd( '/Users/omega/alex/adforesight/' )
source(paste(mypath, '/R/myr2score.R',sep=''))

#output.path <- '/Users/omega/alex/adforesight/mydata/outdata/sd2_projall_noscale/'
if (dir.exists(output.path)){ 1} else {dir.create(output.path, recursive=TRUE)}
#if (dir.exists(paste(output.path,'/',extension, '/', sep='')){ 1} else {dir.create(paste(output.path,'/',extension, '/'), recursive=TRUE)}

# Load in Behavior and Imaging Data
behavior <- read.csv('./mydata/All_Behavior.csv')
behavior<-behavior[c(1:24),];
mygroup <- behavior$genotype[1:24] #superfluous junk
myindex <- c(1:24)


mytemplate<-antsImageRead('/Users/alex/brain_data/MEMRI_CVN_NOS/MEMRI_RARE/median_images/MDT_T2.nii.gz')
mask  <- antsImageRead('./mydata/MDT_mask_e3.nii')
mask <- thresholdImage( mask, 0.1, Inf )

mask_pos_jac <-antsImageRead('/Users/alex/brain_data/MEMRI_CVN_NOS/MEMRI_VBA/jac_ad_gt_ctrl_fdrp0p01k100.nii.gz');
mask_neg_jac<--antsImageRead('/Users/alex/brain_data/MEMRI_CVN_NOS/MEMRI_VBA/jac_ctrl_gt_ad_fdr0p01k100.nii.gz');
mask_neg_jac<-antsImageRead('/Users/alex/brain_data/MEMRI_CVN_NOS/MEMRI_VBA/spm_results/JAC/jac_ctrl_gt_ad_fdr0p01k100.nii')
mask_pos_jac<-antsImageRead('/Users/alex/brain_data/MEMRI_CVN_NOS/MEMRI_VBA/spm_results/JAC/jac_ad_gt_ctrl_fdrp0p01k100.nii')

mask_pos_jac<-thresholdImage( mask_pos_jac, 0.0001, Inf )
mask_neg_jac<-thresholdImage( mask_neg_jac, 0.0001, Inf )

#mask_pos_Mn <-antsImageRead('/Users/alex/brain_data/MEMRI_CVN_NOS/MEMRI_VBA/Mn_ad_gt_ctrl_0p2k100.nii.gz');
#mask_neg_Mn<--antsImageRead('/Users/alex/brain_data/MEMRI_CVN_NOS/MEMRI_VBA/Mn_ctrl_gt_ad_0p2k100.nii.gz');corrupt
##
mask_pos_Mn <-antsImageRead('/Users/alex/brain_data/MEMRI_CVN_NOS/MEMRI_VBA/spm_results/T1/T1_ad_gt_ctrl_0p2k100.nii')
mask_neg_Mn<-antsImageRead('/Users/alex/brain_data/MEMRI_CVN_NOS/MEMRI_VBA/spm_results/T1/T1_ctrl_gt_ad_0p2k100.nii')
##
mask_pos_Mn<-thresholdImage( mask_pos_Mn, 0.0001, Inf )
mask_neg_Mn<-thresholdImage( (mask_neg_Mn), 0.0001, Inf )

mask_pos_chi <-antsImageRead('/Users/alex/brain_data/MEMRI_CVN_NOS/MEMRI_VBA/Chi_ad_gt_ctrl_k100_fdr0p1.nii.gz');
mask_neg_chi<--antsImageRead('/Users/alex/brain_data/MEMRI_CVN_NOS/MEMRI_VBA/Chi_ctrl_gt_ad_k100_fdr0p1.nii.gz');
mask_neg_chi<-antsImageRead('/Users/alex/brain_data/MEMRI_CVN_NOS/MEMRI_VBA/spm_results/Chi_subtracted/CHI_ctrl_gt_ad_k100_0p1.nii')
mask_pos_chi <-antsImageRead('/Users/alex/brain_data/MEMRI_CVN_NOS/MEMRI_VBA/spm_results/Chi_subtracted/CHI_AD_gt_ctrl_k100_fdr0p1.nii')

mask_pos_chi<-thresholdImage( mask_pos_chi, 0.0001, Inf )
mask_neg_chi<-thresholdImage( mask_neg_chi, 0.0001, Inf )
#plot(mytemplate, list(mask_neg_chi))
#plot(mytemplate, list(mask_pos_chi))
###

mang_files <- list.files(path = "./mydata/imdata/", pattern = "T2_to_MDT",full.names = T,recursive = T)
mang_pos_Mn<- imagesToMatrix(mang_files,mask_pos_Mn)
mang_neg_Mn<- imagesToMatrix(mang_files,mask_neg_Mn)
mask_combo_Mn<-imageListToMatrix(list(mask_pos_Mn,mask_neg_Mn),mask = mask)
#sanity check
#plot(mytemplate, list(mask_pos_Mn,mask_neg_Mn))

jac_files <- list.files(path = "./mydata/imdata/", pattern = "jac_to_MDT",full.names = T,recursive = T)
mang_pos_jac<- imagesToMatrix(jac_files,mask_pos_jac)
mang_neg_jac<- imagesToMatrix(jac_files,mask_neg_jac)
mask_combo_jac<-imageListToMatrix(list(mask_pos_jac,mask_neg_jac),mask = mask)
#sanity check
#plot(mytemplate, list(mask_pos_jac,mask_neg_jac))

suscept_files <- list.files(path = "./mydata/imdata/", pattern = "X_to_MDT",full.names = T,recursive = T)
mang_pos_chi<- imagesToMatrix(suscept_files,mask_pos_chi)
mang_neg_chi<- imagesToMatrix(suscept_files,mask_neg_chi)
mask_combo_chi<-imageListToMatrix(list(mask_pos_chi,mask_neg_chi),mask = mask)
#sanity check
#plot(mytemplate, list(mask_pos_chi,mask_neg_chi))

#pick your contrast 
mang_files <- jac_files
mang_mat <- imagesToMatrix(mang_files,mask)

###############################################
###make models
mygroup <- behavior$genotype[1:24]
myindex <- c(1:24)

mydfb <- data.frame("mysubject_index" = factor(as.integer(myindex)),"mygenotype"=mygroup)
kable(mydfb, align = 'c')

set.seed(1)
k<-4

performances <- c()
myBICs  <- c()
myR2score<-c()
myps<-c()
gfit<-c()

###build k models and retain the best performing one in terms of RMSE2
#LOOCV replaces folds
#k<-length(rows.train)-1

k<-4
set.seed(1)
res_train<-createFolds(behavior$genotype,k, list = TRUE, returnTrain = TRUE)
set.seed(1)
res_test<-createFolds(behavior$genotype,k)


for (myfold in 1:k){
  # for (myfold in 3){
  gc(verbose = TRUE, reset = FALSE)
  print('myfold:',myfold)
  print(myfold)
  
  rows.train<-as.integer(unlist(res_train[myfold]))
  rows.test<-as.integer(unlist(res_test[myfold]))
  
  #pos clusters only
  mang_pos.train <- mang_pos_Mn[rows.train, ]
  mang_pos.test <- mang_pos_Mn[rows.test, ]
  imgpredtrain_mang<-mang_pos.train  %*% t(imageListToMatrix(list(mask)))
  
  #neg_clusters only
  mang_neg.train <- mang_neg_Mn[rows.train, ]
  mang_neg.test <- mang_neg_Mn[rows.test, ]
  
  #combo clusters
  
  
  behav.train <- behavior[rows.train, ]
  behav.test  <- behavior[rows.test, ]
  dist4.train <- behav.train[,'d4']
  dist4.test <- behav.test[,'d4'] 
  
}


###end make models
###################################

#back <-matrixToImages(jac_mat,mask)

jac.train <- jac_mat[rows.train, ]
jac.test <- jac_mat [rows.test, ]
behav.train <- behavior[rows.train, ]
behav.test <-behavior[rows.test, ]
imgpredtrain_jac<-jac.train  %*% t(imageListToMatrix(list(mask)))
imgpredtest_jac<-jac.test # %*% t(mask)

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


#Validation part
rows.valid<-c(1:24_)
mang.valid <- mang_mat[rows.valid, ]

#read eigenregions for best myfold
#paste(output.path,extension,'sd2eig' ,as.character(i), 'fold', toString(myfold), '.nii.gz',sep='')

eig_files_Mn <- list.files(path = paste(vbaoutput.path,sep=''), pattern=paste('*', 'Mn', sep=''),full.names = T,recursive = T)
eig_mat_Mn <- imagesToMatrix(eig_files_Mn,mask)
imgmat_mang_valid <- mang.valid %*% t(eig_mat_Mn) # [24,numvox] [nvecsx3,numvox]


dist4.valid  <- behavior[rows.valid, 'd4']

projs.valid <- data.frame(cbind(dist4.valid,imgmat_mang_valid))
colnames(projs.valid) <- c('Dist_4', paste0('Proj', c(1:ncolcombo)))
distpred <- predict.lm(mylm, newdata=projs.valid) 
mymodel<-lm(distpred~dist4.valid)

RSS <- c(crossprod(mymodel$residuals))
MSE <- RSS / length(mymodel$residuals)
RMSE <- sqrt(MSE)
RMSE2<-sqrt(mean((distpred - dist4.valid)^2)) 
mysummary <-summary(mymodel)
r2pred <- mysummary$adj.r.squared
ppred <- mysummary$coefficients[2,4]


max(behavior$d4[1:24])
RMSE_valid<-RMSE2
BIC_valid <- BIC(mymodel)
R2score_valid<-myr2score(distpred,dist4.valid)
res_cor<-cor.test(dist4.valid,distpred)

myplot<- visreg(mymodel, gg=TRUE, scale='linear', plot=TRUE, xlim=c(0,max(dist4.valid)),ylim=c(0,max(dist4.valid))) 
myplot2<-plot(myplot,xlim=c(0,max(dist4.valid)),ylim=c(0,max(dist4.valid)))
ggsave(paste(output.path,extension,'MnValidationSet',toString(myfold),'sd2plainjane.pdf',sep=''), plot = last_plot(), device = 'pdf', 
       scale = 1, width = 4, height = 4, units = c("in"),dpi = 300)

myplot + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
               panel.background = element_rect(fill = "transparent", colour = NA),
               axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
               axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))+
  #xlim(0,1200)+ylim(0,1200)+coord_cartesian(xlim = c(1200, 1200),ylim = c(1200, 1200)) + coord_equal()+
  ggtitle(paste("RMSE=",formatC(RMSE_valid,digits=2, format="f"), 
                # "R2score=",formatC(R2score_valid,digits=2, format="f"), 
                # " R2=", formatC(r2pred,digits=2, format="f"), 
                " p= ", formatC(ppred,digits=4, format="f"),
                "  BIC=", formatC(BIC_valid,digits=2, format="f")))

ggsave(paste(output.path,extension,'MnValidationSet',toString(myfold),'sd2.pdf',sep=''), plot = last_plot(), device = 'pdf', 
       scale = 1, width = 4, height = 4, units = c("in"),dpi = 300)

numcols<-dim(projs.valid)[2]
rd4 <- t(t(dist4.valid)[rep(1,c(3*mynvecs)),])
pcor<-c(numcols)
corval<-c(numcols)


for (i in 1:numcols) {
  mypcor<-cor.test(t(dist4.valid),t(projs.valid[,i]))
  pcor[i]<-mypcor$p.value
  corval[i]<-mypcor$estimate
}

rt<-glance(cor.test(dist4.valid,distpred))
corval[1]<-rt$estimate
pcor[1]<-rt$p.value
# for (i in 1:numcols){
#   res1eig<-matrixToImages(jeanat_mang$fusedlist,mask = mask)[[i]] #eig1
#   antsImageWrite(res1eig,paste(output.path,extension,'JoinEanat' ,as.character(i), '.nii.gz',sep=''))
# }

mycorsdf_eig2d4<-data.frame(rbind(pcor,corval),row.names=c("pcor","cor"))
colnames(mycorsdf_eig2d4)<-c('total', paste0('Proj', c(1:ncolcombo)))
write.csv(mycorsdf_eig2d4, file = paste(output.path ,extension,'fold', toString(myfold), 'd4corsvalidsd2.csv',sep=''))

myperf<-data.frame(rbind(distpred,dist4.valid),row.names=c("d_predicted","d_valid"))
write.csv(myperf, file = paste(output.path ,extension,'fold', toString(myfold), 'distances4_validsd2.csv',sep=''))
