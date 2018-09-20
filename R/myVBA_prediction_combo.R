library( knitr )
library(ANTsR)
library(visreg)
library(robustbase)
library(groupdata2) 
library(ggplot2)
library(caret)
require(broom)

#[18,num vox]  * [numvox,nvecs]
# Load in Behavior and Imaging Data
mypath<-'/Users/alex/GitHub/adforesight/' #flavors of serifos '/Volumes/CivmUsers/omega/alex/GitHub/adforesight/'
# Load in Behavior and Imaging Data
#setwd( '/Users/alex/GitHub/adforesight/' )
mypath<-'/Users/alex/Documents/GitHub/adforesight/' #flavors of ithaka

setwd(mypath)
mysuffix<-'jac'
#mysuffix<-'Mn'
#mysuffix<-'chi'

ncolcombo<-2 #for pos and negative clusters
extension<-paste('VBA_',mysuffix,'_posneg_pred',sep='')
output.path <- paste(mypath,'/mydata/outdata/',extension, '/', sep='') #sd2_projall_noscale/'
vbainput.path<-'/Users/alex/Documents/GitHub/adforesight/mydata/vba/' #flavors of ithaka
vbainput.path<-'/Users/alex/brain_data/MEMRI_CVN_NOS/MEMRI_VBA/' #flavors of serifos
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
mydfb <- data.frame("mysubject_index" = factor(as.integer(myindex)),"mygenotype"=mygroup)
kable(mydfb, align = 'c')


#mytemplate<-antsImageRead('/Users/alex/brain_data/MEMRI_CVN_NOS/MEMRI_RARE/median_images/MDT_T2.nii.gz')
mytemplate<-antsImageRead('./mydata/median_images/MDT_T2.nii.gz')
mask  <- antsImageRead('./mydata/MDT_mask_e3.nii')
mask <- thresholdImage( mask, 0.1, Inf )

#flavors of serifos
# mask_pos_jac <-antsImageRead('/Users/alex/brain_data/MEMRI_CVN_NOS/MEMRI_VBA/jac_ad_gt_ctrl_fdrp0p01k100.nii.gz');
# mask_neg_jac<--antsImageRead('/Users/alex/brain_data/MEMRI_CVN_NOS/MEMRI_VBA/jac_ctrl_gt_ad_fdr0p01k100.nii.gz');
# mask_neg_jac<-antsImageRead('/Users/alex/brain_data/MEMRI_CVN_NOS/MEMRI_VBA/spm_results/JAC/jac_ctrl_gt_ad_fdr0p01k100.nii')
# mask_pos_jac<-antsImageRead('/Users/alex/brain_data/MEMRI_CVN_NOS/MEMRI_VBA/spm_results/JAC/jac_ad_gt_ctrl_fdrp0p01k100.nii')
# mask_pos_chi <-antsImageRead('/Users/alex/brain_data/MEMRI_CVN_NOS/MEMRI_VBA/Chi_ad_gt_ctrl_k100_fdr0p1.nii.gz');
# mask_neg_chi<--antsImageRead('/Users/alex/brain_data/MEMRI_CVN_NOS/MEMRI_VBA/Chi_ctrl_gt_ad_k100_fdr0p1.nii.gz');
# mask_neg_chi<-antsImageRead('/Users/alex/brain_data/MEMRI_CVN_NOS/MEMRI_VBA/spm_results/Chi_subtracted/CHI_ctrl_gt_ad_k100_0p1.nii')
# mask_pos_chi <-antsImageRead('/Users/alex/brain_data/MEMRI_CVN_NOS/MEMRI_VBA/spm_results/Chi_subtracted/CHI_AD_gt_ctrl_k100_fdr0p1.nii')
# mask_pos_Mn <-antsImageRead('/Users/alex/brain_data/MEMRI_CVN_NOS/MEMRI_VBA/spm_results/T1/T1_ad_gt_ctrl_0p2k100.nii')
# mask_neg_Mn<-antsImageRead('/Users/alex/brain_data/MEMRI_CVN_NOS/MEMRI_VBA/spm_results/T1/T1_ctrl_gt_ad_0p2k100.nii')
# #mask_pos_Mn <-antsImageRead('/Users/alex/brain_data/MEMRI_CVN_NOS/MEMRI_VBA/Mn_ad_gt_ctrl_0p2k100.nii.gz');
# #mask_neg_Mn<--antsImageRead('/Users/alex/brain_data/MEMRI_CVN_NOS/MEMRI_VBA/Mn_ctrl_gt_ad_0p2k100.nii.gz');corrupt

#flavors of ithaka
#reading the positive and negative clusters rom vba
mask_pos_jac <- antsImageRead(paste('./mydata/vba/jac_pos.nii', sep = ''))
mask_neg_jac <- antsImageRead(paste('./mydata/vba/jac_neg.nii', sep = ''))
mask_jac<-mask_pos_jac+mask_neg_jac

mask_pos_Mn <- antsImageRead(paste('./mydata/vba/Mn_pos.nii', sep = ''))
mask_neg_Mn <- antsImageRead(paste('./mydata/vba/Mn_neg.nii', sep = ''))
mask_Mn<-mask_pos_Mn+mask_neg_Mn

mask_pos_chi <- antsImageRead(paste('./mydata/vba/chi_pos.nii', sep = ''))
mask_neg_chi <- antsImageRead(paste('./mydata/vba/chi_neg.nii', sep = ''))
mask_chi<-mask_pos_chi+mask_neg_chi


mask_pos_jac<-thresholdImage( mask_pos_jac, 0.0001, Inf )
mask_neg_jac<-thresholdImage( mask_neg_jac, 0.0001, Inf )
mask_jac<-thresholdImage( mask_jac, 0.0001, Inf )
mask_pos_Mn<-thresholdImage( mask_pos_Mn, 0.0001, Inf )
mask_neg_Mn<-thresholdImage( (mask_neg_Mn), 0.0001, Inf )
mask_Mn<-thresholdImage( mask_Mn, 0.0001, Inf )
mask_pos_chi<-thresholdImage( mask_pos_chi, 0.0001, Inf )
mask_neg_chi<-thresholdImage( mask_neg_chi, 0.0001, Inf )
mask_chi<-thresholdImage( mask_chi, 0.0001, Inf )
#plot(mytemplate, list(mask_neg_chi))
#plot(mytemplate, list(mask_pos_chi))
###
antsImageWrite(mask_pos_jac,paste('./mydata/vba/mask_pos_jac.nii', sep = ''))
antsImageWrite(mask_neg_jac,paste('./mydata/vba/mask_neg_jac.nii', sep = ''))
antsImageWrite(mask_jac,paste('./mydata/vba/mask_jac.nii', sep = ''))
antsImageWrite(mask_pos_Mn,paste('./mydata/vba/mask_pos_Mn.nii', sep = ''))
antsImageWrite(mask_neg_Mn,paste('./mydata/vba/mask_neg_Mn.nii', sep = ''))
antsImageWrite(mask_Mn,paste('./mydata/vba/mask_Mn.nii', sep = ''))
antsImageWrite(mask_pos_chi,paste('./mydata/vba/mask_pos_chi.nii', sep = ''))
antsImageWrite(mask_neg_chi,paste('./mydata/vba/mask_neg_chi.nii', sep = ''))
antsImageWrite(mask_chi,paste('./mydata/vba/mask_chi.nii', sep = ''))

mang_files <- list.files(path = "./mydata/imdata/", pattern = "T2_to_MDT",full.names = T,recursive = T)
pos_Mn<- imagesToMatrix(mang_files,mask_pos_Mn)
neg_Mn<- imagesToMatrix(mang_files,mask_neg_Mn)
mask_combo_Mn<-imageListToMatrix(list(mask_pos_Mn,mask_neg_Mn),mask = mask)
#paranoid check
all_Mn<-imagesToMatrix(mang_files,mask_Mn)
#sanity check
#plot(mytemplate, list(mask_pos_Mn,mask_neg_Mn))

jac_files <- list.files(path = "./mydata/imdata/", pattern = "jac_to_MDT",full.names = T,recursive = T)
pos_jac<- imagesToMatrix(jac_files,mask_pos_jac)
neg_jac<- imagesToMatrix(jac_files,mask_neg_jac)
mask_combo_jac<-imageListToMatrix(list(mask_pos_jac,mask_neg_jac),mask = mask)
all_jac<-imagesToMatrix(jac_files,mask_jac)
#sanity check
#plot(mytemplate, list(mask_pos_jac,mask_neg_jac))

chi_files <- list.files(path = "./mydata/imdata/", pattern = "X_to_MDT",full.names = T,recursive = T)
pos_chi<- imagesToMatrix(chi_files,mask_pos_chi)
neg_chi<- imagesToMatrix(chi_files,mask_neg_chi)
mask_combo_chi<-imageListToMatrix(list(mask_pos_chi,mask_neg_chi),mask = mask)
all_chi<-imagesToMatrix(chi_files,mask_chi)
#sanity check
#plot(mytemplate, list(mask_pos_chi,mask_neg_chi))


###############################################
###initial values for params of our models

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

#pick your contrast 
#alex gotta flex and change this
if (mysuffix=='jac') {
  myfiles <- jac_files
  mymask_pos<-mask_pos_jac
  mymask_neg<-mask_neg_jac
  mymat_pos <- imagesToMatrix(myfiles,mymask_pos)
  mymat_neg <- imagesToMatrix(myfiles,mymask_neg)
  myregion_pos<-pos_jac
  myregion_neg<-neg_jac
  mymask2_pos<-as.matrix(colMeans(myregion_pos))/colMeans(myregion_pos)
  mymask2_neg<-as.matrix(colMeans(myregion_neg))/colMeans(myregion_neg)
} else if (mysuffix=='Mn'){
  myfiles <- mang_files
  mymask_pos<-mask_pos_Mn
  mymask_neg<-mask_neg_Mn
  mymat_pos <- imagesToMatrix(myfiles,mymask_pos)
  mymat_neg <- imagesToMatrix(myfiles,mymask_neg)
  myregion_pos<-pos_Mn
  myregion_neg<-neg_Mn
  mymask2_pos<-as.matrix(colMeans(myregion_pos))/colMeans(myregion_pos)
  mymask2_neg<-as.matrix(colMeans(myregion_neg))/colMeans(myregion_neg)
}else if (mysuffix=='chi'){
  myfiles <- chi_files
  mymask_pos<-mask_pos_chi
  mymask_neg<-mask_neg_chi
  mymat_pos <- imagesToMatrix(myfiles,mymask_pos)
  mymat_neg <- imagesToMatrix(myfiles,mymask_neg)
  myregion_pos<-pos_chi
  myregion_neg<-neg_chi
  mymask2_pos<-as.matrix(colMeans(myregion_pos))/colMeans(myregion_pos)
  mymask2_neg<-as.matrix(colMeans(myregion_neg))/colMeans(myregion_neg)
}else {print('all hell break loose')}

for (myfold in 1:k){
  # for (myfold in 3){
  gc(verbose = TRUE, reset = FALSE)
  print('myfold:',myfold)
  print(myfold)
  
  rows.train<-as.integer(unlist(res_train[myfold]))
  rows.test<-as.integer(unlist(res_test[myfold]))
  
  #pos clusters only
  mang_pos.train <- mymat_pos[rows.train, ]
  mang_neg.train <- mymat_neg[rows.train, ]
  mang_pos.test <- mymat_pos[rows.test, ]
  mang_neg.test <- mymat_neg[rows.test, ]
  
  imgpredtrain_mang_pos<-mang_pos.train  %*% mymask2_pos #gives NA
  imgpredtrain_mang_neg<-mang_neg.train  %*% mymask2_neg #gives NA
  
  imgpredtest_mang_pos<-mang_pos.test  %*% mymask2_pos #gives NA
  imgpredtest_mang_neg<-mang_neg.test  %*% mymask2_neg #gives NA
  
  behav.train <- behavior[rows.train, ]
  behav.test  <- behavior[rows.test, ]
  dist4.train <- behav.train[,'d4']
  dist4.test <- behav.test[,'d4'] 
  
  
  projs.train <- data.frame(cbind(dist4.train,imgpredtrain_mang_pos,imgpredtrain_mang_neg)) # column combine the behavior wth the projections
  colnames(projs.train) <- c('Dist_4', paste0('Proj', c(1:2*ncol(imgpredtrain_mang_pos)))) # insert column names
  projs.test <- data.frame(cbind(dist4.test,imgpredtest_mang_pos,imgpredtest_mang_neg)) # column combine the behavior wth the projections
  colnames(projs.test) <- c('Dist_4', paste0('Proj', c(1:2*ncol(imgpredtest_mang_pos)))) # insert column names
  
  
  mylm <- lm('Dist_4 ~ .', data=projs.train) # behavior correlation with the number of projections
  distpred <- predict.lm(mylm, newdata=projs.test) # based on the linear model predict the distances for the same day
  modsum <-summary(lm(distpred~dist4.test))
  
  
  mylm <- lm('Dist_4 ~ .', data=projs.train) # behavior correlation with projections
  summylm<-summary(mylm)
  summanovalm<-anova(mylm)
  
  rSquared <- summary(mylm)$r.squared
  pVal <- anova(mylm)$'Pr(>F)'[1]
  
  mylmsummary<-glance(mylm)
  pval1<-mylmsummary$p.value
  
  distpred4 <- predict.lm(mylm, newdata=projs.test) # based on the linear model predict the distances for the same day
  glance(cor.test(projs.test$Dist_4,(distpred4)))
  
  glance(cor.test(distpred4,dist4.test))
  #remove next lines for LOOCV
  mymodel<-lm(distpred4~dist4.test)
  modsum <-summary(mymodel)
  r2 <- modsum$r.squared #modsum$adj.r.squared
  # my.p <- modsum$coefficients[2,4]
  
  RMSE2<-sqrt(mean((distpred4 - dist4.test)^2)) 
  performances[myfold]<-RMSE2
  myR2score[myfold]<-myr2score(distpred4,dist4.test)
  myps[myfold]<-pval1<-mylmsummary$p.value #my.p
  myBICs[myfold] <- BIC(mylm)
  
  
  ###
  mytheme <- theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
                   panel.background = element_rect(fill = "white"))
  myplot<- visreg(mymodel, gg=TRUE)
  myplot + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
                 panel.background = element_rect(fill = "transparent", colour = NA),
                 #xaxs="i", yaxs="i",
                 axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                 axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) +
    ggtitle(paste("RMSE=",formatC(RMSE2,digits=2, format="f"), "p=",formatC(myps[myfold],digits=4, format="f"), " BIC=", formatC(BIC(mymodel),digits=2, format="f")))
  
  ggsave(paste(output.path,extension,'fold',toString(myfold),'.pdf',sep=''), plot = last_plot(), device = 'pdf',
         scale = 1, width = 4, height = 4, units = c("in"),dpi = 300)
  
  save(mylm, file=paste(output.path , "model2", toString(myfold), ".Rdata", sep=''))
  save(mymodel, file=paste(output.path , "behavmodelsd2", toString(myfold), ".Rdata", sep=''))
  myperf<-data.frame(rbind(distpred4,dist4.test),row.names=c("d_predicted","d_valid"))
  write.csv(myperf, file = paste(output.path ,extension,'distances4_pv_fold' , toString(myfold), '.csv',sep=''))
  myperf<-data.frame(c(RMSE2,myR2score[myfold],myps[myfold],myBICs[myfold], r2),row.names=c("RMSE2","R2score","p","BIC", "R2"))
  write.csv(myperf, file = paste(output.path ,extension,'distances4_stats_fold' , toString(myfold), '.csv',sep=''))
  
  #plot(dist4.test, distpred4, xlab = 'Real Swim Distance on Day 4', ylab = 'Predicted Swim Distance on day 4', 
  #      main='Predicted vs. Real Swim Distance on Day 4', ylim = c(0,1200),xlim = c(0,1200)) # generate plot
  # 
  
  
}


###end make models
###################################

###################################
#####    for validation now    ####
###################################

myminfold<-which(performances == min(performances), arr.ind = TRUE)
myfold<-myminfold
load(file=paste(output.path , "model2", toString(myfold), ".Rdata", sep='')) # loads mylm

rows.valid <- c(1:24)


mang.valid_pos <- mymat_pos[rows.valid, ]
mang.valid_neg <- mymat_neg[rows.valid, ]
imgmat_mang_valid_pos <- mang.valid_pos %*% mymask2_pos
imgmat_mang_valid_neg <- mang.valid_neg %*% mymask2_neg

dist4.valid  <- behavior[rows.valid, 'd4']

projs.valid <- data.frame(cbind(dist4.valid,imgmat_mang_valid_pos,imgmat_mang_valid_neg))
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
rd4 <- t(t(dist4.valid)[rep(1,c(numcols)),])
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

