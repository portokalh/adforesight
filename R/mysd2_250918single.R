#load libs
library( knitr )
library(ANTsR)
library(visreg)
library(robustbase)
library(groupdata2) 
library(ggplot2)
library(caret)
require(broom)

#Manganese Enhanced MRI   predicts cognitive  performnace Alex Badea and Natalie Delpratt
#reuses results of sparsedecom2 for best performing fold Alex Badea 7 dept 2018
#uses RMSE as prediction error, rather than goodness of fit error 30 August 2018

#legacy from Natalie to rememebr paths
#setwd( '/Users/omega/alex/adforesight/' )
#output.path <- '/Users/omega/alex/adforesight/mydata/outdata/sd2_projall_noscale/' #ND


mypath<-'/Volumes/CivmUsers/omega/alex/GitHub/adforesight/'
mypath <- '/Users/alex/GitHub/adforesight/' #flavors of serifos
mypath <- '/Users/alex/Documents/GitHub/adforesight/' #flavors of ithaka
setwd(mypath)
source(paste(mypath, '/R/myr2score.R',sep=''))

mysp <- 0.05 #0.05  # 0.01  # 0.05 #0.2 #0.05 # was 0.005 sparseness
#if mynvecs is set to one eig1 neetds to be transposed
mynvecs <- 2 # 5 vecs is better # 10 shows ventral thalamic nuclei 10 # 50 #put back to 50 alex #tested with 10

myell1 <- 1 # make this smaller 0.5, Brian says it is not what i think just switch between l0 and l1
myits<-15 #15 #5 
mysmooth<-0.1 #0.1 # 0.1 #0.01 #0  # was 0.01
myclus<-250 #was 250

# Load in Behavior and Imaging Data
behavior <- read.csv('./mydata/All_Behavior.csv')
labled.set <-read.csv('./mydata/legendsCHASS2symmetric.csv')
labeled.brain.img <-  antsImageRead('./mydata/MDT_labels_chass_symmetric.nii.gz')
mask  <- antsImageRead('./mydata/MDT_mask_e3.nii')
mask <- thresholdImage( mask, 0.1, Inf )

#read all 3 contrast files
mang_files <- list.files(path = "./mydata/imdata/", pattern = "T2_to_MDT",full.names = T,recursive = T)
jac_files <- list.files(path = "./mydata/imdata/", pattern = "jac_to_MDT",full.names = T,recursive = T)
chi_files <- list.files(path = "./mydata/imdata/", pattern = "X_to_MDT",full.names = T,recursive = T)


########################################
#build a place to save results
extension<-paste('sd2SINGLEJAC', 'sp', toString(mysp), 'vecs', toString(mynvecs), 's', toString(mysmooth),'clus', toString(myclus), sep='') # 'JACsp0p005s0'
output.path <- paste(mypath,'/mydata/outdata_sd2/',extension, '/', sep='') #sd2_projall_noscale/'
if (dir.exists(output.path)){ 1} else {dir.create(output.path, recursive=TRUE)}

#pick yourcontrast
mang_mat <- imagesToMatrix(jac_files,mask)
#######################################
#let things flow from here

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
#considet using LOOCV to replace folds, but results may be unstable
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
  
  mang.train <- mang_mat[rows.train, ]
  mang.test <- mang_mat[rows.test, ]
  
  
  behav.train <- behavior[rows.train, ]
  behav.test  <- behavior[rows.test, ]
  dist4.train <- behav.train[,'d4']
  dist4.test <- behav.test[,'d4'] 
  
  start_time <- Sys.time()
 #negative sparseness is what? allows for negative weights!
  myeig2_mang<-sparseDecom2(inmatrix = list(mang.train,as.matrix(behav.train$d4)),its = myits, cthresh=c(myclus,0), smooth = mysmooth, mycoption = 0, sparseness = c(mysp,1), nvecs = mynvecs, verbose=1, statdir=paste(output.path))
  
  
  #myeig2_mang<-sparseDecom(inmatrix = mang.train,its = myits, cthresh=c(myclus), smooth = mysmooth, mycoption = 0, sparseness = c(mysp), nvecs = mynvecs, verbose=1, statdir=paste(output.path2))
   
  end_time <- Sys.time()
  t1time<-end_time - start_time
  print(t1time)
  
  imgpredtrain_mang<-mang.train %*% (myeig2_mang$eig1)
  imgpredtest_mang<-mang.test %*% (myeig2_mang$eig1)
  
  
  ####start do single alex
  ncolcombo<-ncol( imgpredtrain_mang)
  
  projs.train <- data.frame(dist4.train, imgpredtrain_mang) # column combind the behavior wth the projections
  colnames(projs.train) <- c('Dist_4', paste0('Proj', c(1:ncolcombo)))
  projs.test <- data.frame(dist4.test, imgpredtest_mang ) # column combind the behavior wth the projections
  colnames(projs.test) <- c('Dist_4', paste0('Proj', c(1:ncolcombo)))
  ###end do single alex
  
  
  mylm <- lm('Dist_4 ~ .', data=projs.train) # behavior correlation with projections
  summylm<-summary(mylm)
  summanovalm<-anova(mylm)
  
  rSquared <- summary(mylm)$r.squared
  pVal <- anova(mylm)$'Pr(>F)'[1]

  mylmsummary<-glance(mylm)
  pval1<-mylmsummary$p.value
  
  e2i_mang<-matrixToImages((t(myeig2_mang$eig1)),mask = mask)
   
  
  for (i in 1:mynvecs){
    antsImageWrite(e2i_mang[[i]],paste(output.path,extension,'sd2eig_' ,as.character(i), 'fold', toString(myfold), '_Mn.nii.gz',sep=''))
    # antsImageWrite(e2i_jac[[i]],paste(output.path,extension,'sd2eig_' ,as.character(i), 'fold', toString(myfold), '_jac.nii.gz',sep=''))
    # antsImageWrite(e2i_chi[[i]],paste(output.path,extension,'sd2eig_' ,as.character(i), 'fold', toString(myfold), '_chi.nii.gz',sep=''))
    }
  
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

  ggsave(paste(output.path,extension,'Mnfold',toString(myfold),'.pdf',sep=''), plot = last_plot(), device = 'pdf',
         scale = 1, width = 4, height = 4, units = c("in"),dpi = 300)

  save(mylm, file=paste(output.path , "model2", toString(myfold), ".Rdata", sep=''))
  save(mymodel, file=paste(output.path , "behavmodelsd2", toString(myfold), ".Rdata", sep=''))
  myperf<-data.frame(rbind(distpred4,dist4.test),row.names=c("d_predicted","d_valid"))
  write.csv(myperf, file = paste(output.path ,extension,'distances4_pv_fold' , toString(myfold), '.csv',sep=''))
  myperf<-data.frame(c(RMSE2,myR2score[myfold],myps[myfold],myBICs[myfold], r2),row.names=c("RMSE2","R2score","p","BIC", "R2"))
  write.csv(myperf, file = paste(output.path ,extension,'distances4_stats_fold' , toString(myfold), '.csv',sep=''))
  gc(verbose = TRUE, reset = FALSE)
 
  }


###################################
#####    for validation now    ####
###################################

myminfold<-which(performances == min(performances), arr.ind = TRUE)
myfold<-myminfold
load(file=paste(output.path , "model2", toString(myfold), ".Rdata", sep='')) # loads mylm
ncolcombo<-mynvecs
rows.valid <- c(1:24)


mang.valid <- mang_mat[rows.valid, ]

#read eigenregions for best myfold
#paste(output.path,extension,'sd2eig' ,as.character(i), 'fold', toString(myfold), '.nii.gz',sep='')

eig_files_Mn <- list.files(path = paste(output.path,sep=''), pattern=paste('*', 'fold', toString(myfold), '_Mn.nii.gz', sep=''),full.names = T,recursive = T)
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


mycorsdf_eig2d4<-data.frame(rbind(pcor,corval),row.names=c("pcor","cor"))
colnames(mycorsdf_eig2d4)<-c('total', paste0('Proj', c(1:ncolcombo)))
write.csv(mycorsdf_eig2d4, file = paste(output.path ,extension,'fold', toString(myfold), 'd4corsvalidsd2.csv',sep=''))

myperf<-data.frame(rbind(distpred,dist4.valid),row.names=c("d_predicted","d_valid"))
write.csv(myperf, file = paste(output.path ,extension,'fold', toString(myfold), 'distances4_validsd2.csv',sep=''))

myeig2_mang_valid<-sparseDecom2(inmatrix = list(mang.valid,as.matrix(dist4.valid)),its = myits, cthresh=c(myclus,0), smooth = mysmooth, mycoption = 0, sparseness = c(mysp,1), nvecs = mynvecs, verbose=1, statdir=paste(output.path))
#myeig2_mang<-sparseDecom2(inmatrix = list(mang.train,as.matrix(behav.train$d4)),its = myits, cthresh=c(myclus,0), smooth = mysmooth, mycoption = 0, sparseness = c(mysp,1), nvecs = mynvecs, verbose=1, statdir=paste(output.path))

e2i_mang_valid<-matrixToImages((t(myeig2_mang_valid$eig1)),mask = mask)

#imgpredtrain_mang<-mang.train %*% (myeig2_mang$eig1)

for (i in 1:mynvecs){
  antsImageWrite(e2i_mang_valid[[i]],paste(output.path,extension,'full_eig' ,as.character(i), '_Mn.nii.gz',sep=''))
  }

#redo fold min or save models

imgmat_mang_valid <- mang.valid %*% t(e2i_mang_valid) # [24,numvox] [nvecsx3,numvox]
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
ggsave(paste(output.path,extension,'newMnValidationSetFULL',toString(myfold),'sd2plainjane.pdf',sep=''), plot = last_plot(), device = 'pdf', 
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

ggsave(paste(output.path,extension,'newMnValidationSetFULL',toString(myfold),'sd2.pdf',sep=''), plot = last_plot(), device = 'pdf', 
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


mycorsdf_eig2d4<-data.frame(rbind(pcor,corval),row.names=c("pcor","cor"))
colnames(mycorsdf_eig2d4)<-c('total', paste0('Proj', c(1:ncolcombo)))
write.csv(mycorsdf_eig2d4, file = paste(output.path ,extension,'FULL', 'd4corsvalidsd2.csv',sep=''))

myperf<-data.frame(rbind(distpred,dist4.valid),row.names=c("d_predicted","d_valid"))
write.csv(myperf, file = paste(output.path ,extension,'FULL', 'distances4_validsd2.csv',sep=''))
