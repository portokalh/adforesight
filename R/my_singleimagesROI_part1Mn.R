##Prior Based Prediction using combined imaging modalities & SparseDecom2

library( knitr )
library(ANTsR)
library(visreg)
library(robustbase)
library(groupdata2) 
library(ggplot2)
library(caret)
require(broom)

mypath<-'/Users/alex/GitHub/adforesight/' #flavors of serifos '/Volumes/CivmUsers/omega/alex/GitHub/adforesight/'
mypath<-'/Users/alex/Documents/GitHub/adforesight/' #flavors of ithaka
setwd(mypath)

Mn_files <- list.files(path = "./mydata/imdata/", pattern = "T2_to_MDT",full.names = T,recursive = T)
jac_files <- list.files(path = "./mydata/imdata/", pattern = "jac_to_MDT",full.names = T,recursive = T)
chi_files <- list.files(path = "./mydata/imdata/", pattern = "X_to_MDT",full.names = T,recursive = T)
Mn_mat <-imagesToMatrix(Mn_files,mask)
jac_mat <-imagesToMatrix(jac_files,mask)
chi_mat <-imagesToMatrix(chi_files,mask)

mysuffix<-'single_ROI'
mycontrast<-'Mn'
my_mat<-Mn_mat
my_files<-Mn_files

mysp <- 0.05  # 0.01  # 0.05 #0.2 #0.05 # was 0.005 sparseness
#if mynvecs is set to one eig1 neetds to be transposed
mynvecs <- 1 # 5 vecs is better # 10 shows ventral thalamic nuclei 10 # 50 #put back to 50 alex #tested with 10
myell1 <- 1 # just switch between l0 (if neg) and l1 (if pos)
myits<-15 # 15
mysmooth<-0.1 #0.1 # 0.1 #0.01 #0  # was 0.01
myclus<-250 

ncolcombo<-mynvecs #mutliply by 3 if all contrasts 
extension<-paste(mysuffix,'_pred',sep='')
output.path0 <- paste(mypath,'/mydata/outdata/',extension, 'sp', toString(mysp),'smooth', toString(mysmooth), '/', sep='') #sd2_projall_noscale/'
if (dir.exists(output.path0)){1} else {dir.create(output.path0, recursive=TRUE)}
source(paste(mypath, '/R/myr2score.R',sep=''))

# Load in Behavior and Imaging Data
behavior <- read.csv('./mydata/All_Behavior.csv')
behavior<-behavior[c(1:24),];
mygroup <- behavior$genotype[1:24] #superfluous junk
myindex <- c(1:24)
mydfb <- data.frame("mysubject_index" = factor(as.integer(myindex)),"mygenotype"=mygroup)
kable(mydfb, align = 'c')



labeled.set <-read.csv('./mydata/legendsCHASS2symmetric.csv')
labeled.brain.img <-  antsImageRead('./mydata/MDT_labels_chass_symmetric.nii.gz')
mytemplate<-antsImageRead('./mydata/median_images/MDT_T2.nii.gz')
mask  <- antsImageRead('./mydata/MDT_mask_e3.nii')
mask <- thresholdImage( mask, 0.1, Inf )





#######################################
#let things flow from here

mygroup <- behavior$genotype[1:24]
myindex <- c(1:24)

mydfb <- data.frame("mysubject_index" = factor(as.integer(myindex)),"mygenotype"=mygroup)
kable(mydfb, align = 'c')

set.seed(1)
k<-4 # 4 fold CV

performances <- c()
myBICs  <- c()
myR2score<-c()
myps<-c()
gfit<-c()


#pick a few regions of interest with anatomical labels
label.ind <- c(118, 20, 19,42,47,51,54,57,59,64,65,120, 121,122, 123) #20,62,63,73,81,83,91,91,119; c(62,122, 67) #

#Both (Left and Right) ROIS
#for ( j in 2:length(label.ind)){
 # for ( j in 1:length(label.ind)){
    for ( j in 1:1){
  mylabel<-as.character(labeled.set$Abbreviation[label.ind[j]])
  #overwrite outputh path now
  output.path <- paste(output.path0, '/', mycontrast, '/', mylabel, '/', sep='') #sd2_projall_noscale/'
  
  if (dir.exists(output.path)){ 1} else {dir.create(output.path, recursive=TRUE)}
  labels<-c(label.ind[j],label.ind[j]+1000)
  initmat<-matrix( rep(0,sum(mask==1)*length(labels)), nrow=length(labels) )
  output.file <- paste(output.path, as.character(labeled.set$Abbreviation[label.ind[j]]),'LR_all.png',sep = "_")
  output.file3D<-paste(output.path, as.character(labeled.set$Abbreviation[label.ind[j]]),'LR_3Dinit.png',sep = "_")
  
  #fill the matrix with the aal region locations
  for ( i in 1:length(labels) )
  {
    vec<-( labeled.brain.img[ mask == 1 ] == labels[i] )
    initmat[i,]<-as.numeric( vec )
  }
 
  ccainit<-initializeEigenanatomy(initmat,mask)
  
  #brain <- renderSurfaceFunction( surfimg = list( labeled.brain.img ), alphasurf = 0.1,
  #                                funcimg = ccainit$initlist, smoothsval = 0.5, smoothfval = 0.2,
  #                                mycol = c('red','blue') )
  mask2<-ccainit$initlist[[1]]+ccainit$initlist[[2]]
 
  
  k<-4
  set.seed(1)
  res_train<-createFolds(behavior$genotype,k, list = TRUE, returnTrain = TRUE)
  set.seed(1)
  res_test<-createFolds(behavior$genotype,k)
  my_mat <-imagesToMatrix(my_files,mask2)
  
#####start alex folds
  for (myfold in 1:k){
    # for (myfold in 3){
    gc(verbose = TRUE, reset = FALSE)
    print('myfold:',myfold)
    print(myfold)
    
    rows.train<-as.integer(unlist(res_train[myfold]))
    rows.test<-as.integer(unlist(res_test[myfold]))
    
    behav.train <- behavior[rows.train, ]
    behav.test  <- behavior[rows.test, ]
    dist4.train <- behav.train[,'d4']
    dist4.test <- behav.test[,'d4'] 
    
    my.train <- my_mat[rows.train, ]
    my.test <- my_mat[rows.test, ]
    
    
    
    start_time <- Sys.time()
    #negative sparseness is what? allows for negative weights!
    # myeig2_jac<-sparseDecom2(inmatrix = list(jac.train,as.matrix(behav.train$d4)),inmask=list(mask2 ,NA), its = myits, cthresh=c(myclus,0), smooth = mysmooth, mycoption = 0, sparseness = c( mysp, -1), nvecs = mynvecs, verbose=1, statdir=paste(output.path),initializationList = ccainit$initlist, priorWeight = 0.9)
    myeig2<-sparseDecom2(inmatrix = list(my.train,as.matrix(behav.train$d4)),inmask=list(mask2 ,NA), its = myits, cthresh=c(myclus,0), smooth = mysmooth, mycoption = 0, sparseness = c( mysp, -1), nvecs = mynvecs, verbose=1, statdir=paste(output.path),initializationList = list(mask2), priorWeight = 0.9)
    ###blurbs start here
    e2i<-matrixToImages(t((myeig2$eig1)),mask = mask2)
    
    
    # for (i in 1:mynvecs){
    #   #antsImageWrite(ccainit$initlist[[i]],paste(output.path,extension,'ccainit_' ,as.character(i), 'fold', toString(myfold), '_jac.nii.gz',sep=''))
    #   antsImageWrite((e2i_jac)[[i]],paste(output.path,extension,'sd2eig_' ,as.character(i), 'fold', toString(myfold), '_jac.nii.gz',sep=''))
    #   antsImageWrite((e2i_Mn)[[i]],paste(output.path,extension,'sd2eig_' ,as.character(i), 'fold', toString(myfold), '_Mn.nii.gz',sep=''))
    #   antsImageWrite((e2i_chi)[[i]],paste(output.path,extension,'sd2eig_' ,as.character(i), 'fold', toString(myfold), '_chi.nii.gz',sep=''))
    #   
    #   }
    #  
      
    end_time <- Sys.time()
    t1time<-end_time - start_time
    print(t1time)
    
    ###blurbs end here
    
    imgpredtrain<-my.train %*% (myeig2$eig1)
    imgpredtest<-my.test %*% (myeig2$eig1)
    
    ####start do combo alex
    #ncolcombo<-nncol( imgpredtrain_jac)
    projs.train <- data.frame(cbind(dist4.train, imgpredtrain )) # column combind the behavior wth the projections
    colnames(projs.train) <- c('Dist_4', paste0('Proj', c(1:ncolcombo)))
    projs.test <- data.frame(cbind(dist4.test, imgpredtest)) # column combind the behavior wth the projections
    colnames(projs.test) <- c('Dist_4', paste0('Proj', c(1:ncolcombo)))
    
    ###end do combo alex
    
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
  
  
  ####end alex folds
  
  ###################################
  #####    for validation now    ####
  ###################################
  
  myminfold<-which(performances == min(performances), arr.ind = TRUE)
  myfold<-myminfold
  load(file=paste(output.path , "model2", toString(myfold), ".Rdata", sep='')) # loads mylm
  #ncolcombo<-mynvecs
  rows.valid <- c(1:24)
  behavior <- read.csv('./mydata/All_Behavior.csv')
  behavior<-behavior[c(1:24),];
  # my_files <- list.files(path = "./mydata/imdata/", pattern = "jac_to_MDT",full.names = T,recursive = T)
  # my_mat <-imagesToMatrix(my_files,mask2)
  my.valid <- my_mat[rows.valid, ]
  
  dist4.valid  <- behavior[rows.valid, 'd4']
  
  myeig2<-sparseDecom2(inmatrix = list(my.valid,as.matrix(dist4.valid)),inmask=list(mask2 ,NA), its = myits, cthresh=c(myclus,0), smooth = mysmooth, mycoption = 0, sparseness = c( mysp, -1), nvecs = mynvecs, verbose=1, statdir=paste(output.path),initializationList = list(mask2), priorWeight = 0.9)
   
  ###blurbs start here
  
  # e2i_Mn<-matrixToImages(t((myeig2_Mn$eig1)),mask = mask2)
  # e2i_jac<-matrixToImages(t((myeig2_jac$eig1)),mask = mask2)
  # e2i_chi<-matrixToImages(t((myeig2_chi$eig1)),mask = mask2)
  # 
  
  
    imgmat_valid <- my.valid %*% (myeig2$eig1)
  
  
  
  projs.valid <- data.frame(cbind(dist4.valid,imgmat_valid))
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
  
} #end for all j labels