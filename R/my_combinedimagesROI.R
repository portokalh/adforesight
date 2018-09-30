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

mysuffix<-'combo_ROI'
mysp <- 0.05 #0.05  # 0.01  # 0.05 #0.2 #0.05 # was 0.005 sparseness
#if mynvecs is set to one eig1 neetds to be transposed
mynvecs <- 1 # 5 vecs is better # 10 shows ventral thalamic nuclei 10 # 50 #put back to 50 alex #tested with 10
myell1 <- 1 # just switch between l0 (if neg) and l1 (if pos)
myits<-15 # 
mysmooth<-0.1 #0.1 # 0.1 #0.01 #0  # was 0.01
myclus<-250 

ncolcombo<-3 #for one region we have 3 contrasts
extension<-paste(mysuffix,'_pred',sep='')
output.path <- paste(mypath,'/mydata/outdata/',extension, 'sp', toString(mysp),'smooth', toString(mysmooth), '/', sep='') #sd2_projall_noscale/'
if (dir.exists(output.path)){1} else {dir.create(output.path, recursive=TRUE)}
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

Mn_files <- list.files(path = "./mydata/imdata/", pattern = "T2_to_MDT",full.names = T,recursive = T)
jac_files <- list.files(path = "./mydata/imdata/", pattern = "jac_to_MDT",full.names = T,recursive = T)
chi_files <- list.files(path = "./mydata/imdata/", pattern = "X_to_MDT",full.names = T,recursive = T)
Mn_mat <-imagesToMatrix(Mn_files,mask)
jac_mat <-imagesToMatrix(jac_files,mask)
chi_mat <-imagesToMatrix(chi_files,mask)


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
label.ind <- c(19,20,42,47,51,54,57,59,62,63,64,65,73,81,83,91,119,121,123)

#Both (Left and Right) ROIS
#for ( j in 2:length(label.ind)){
  for ( j in 5:6){
  mylabel<-as.character(labeled.set$Abbreviation[label.ind[j]])
  #overwrite outputh path now
  output.path <- paste(mypath,'/mydata/outdata/',extension, '/', mylabel, '/', sep='') #sd2_projall_noscale/'
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
  # newmat <- as.matrix(colSums(initmat))
  # newmat <- t(newmat)
#   mask2<-mask[labeled.brain.img[ mask == 1 ] == labels[1]] 
#   mask2<-labeled.brain.img == labels[1]]
# res1<-mask[labeled.brain.img==labels[1]]
#   
#   antsImageWrite(mask2,paste(output.path,extension,'ccainit_mask' ,as.character(i), 'fold', toString(myfold), '_Mn.nii.gz',sep=''))
#   
  ccainit<-initializeEigenanatomy(initmat,mask)
  
  brain <- renderSurfaceFunction( surfimg = list( labeled.brain.img ), alphasurf = 0.1,
                                  funcimg = ccainit$initlist, smoothsval = 0.5, smoothfval = 0.2,
                                  mycol = c('red','blue') )
  mask2<-ccainit$initlist[[1]]+ccainit$initlist[[2]]
  brai2n <- renderSurfaceFunction( surfimg = list( mask2 ), alphasurf = 0.1, funcimg = ccainit$initlist, smoothsval = 0.5, smoothfval = 0.2, mycol = c('red','blue') )
  
  id <- par3d( "userMatrix" )
  rid <- rotate3d( id, pi / 2, 1, 0, 0 )
  rid2 <- rotate3d( id, pi / 2, 0, 0, 1 )
  par3d( userMatrix = id )
  
  dd <- make3ViewPNG( rid, id, rid2, output.file3D)
  par3d( userMatrix = id )
  
  k<-4
  set.seed(1)
  res_train<-createFolds(behavior$genotype,k, list = TRUE, returnTrain = TRUE)
  set.seed(1)
  res_test<-createFolds(behavior$genotype,k)
  
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
    
    Mn.train <- Mn_mat[rows.train, ]
    Mn.test <- Mn_mat[rows.test, ]
    
    jac.train <- jac_mat[rows.train, ]
    jac.test <- jac_mat[rows.test, ]
    
    chi.train <- chi_mat[rows.train, ]
    chi.test <- chi_mat[rows.test, ]
    
    start_time <- Sys.time()
    #negative sparseness is what? allows for negative weights!
    myeig2_jac<-sparseDecom2(inmatrix = list(jac.train,as.matrix(behav.train$d4)),inmask=list(mask,NA), its = myits, cthresh=c(myclus,0), smooth = mysmooth, mycoption = 0, sparseness = c( mysp, -0.5 ), nvecs = mynvecs, verbose=1, statdir=paste(output.path),initializationList = ccainit$initlist, priorWeight = 0.9)
    myeig2_Mn<-sparseDecom2(inmatrix = list(Mn.train,as.matrix(behav.train$d4)),inmask=list(mask,NA), its = myits, cthresh=c(myclus,0), smooth = mysmooth, mycoption = 0, sparseness = c( mysp, -0.5 ), nvecs = mynvecs, verbose=1, statdir=paste(output.path),initializationList = ccainit$initlist, priorWeight = 0.9)
    myeig2_chi<-sparseDecom2(inmatrix = list(chi.train,as.matrix(behav.train$d4)),inmask=list(mask,NA), its = myits, cthresh=c(myclus,0), smooth = mysmooth, mycoption = 0, sparseness = c( mysp, -0.5 ), nvecs = mynvecs, verbose=1, statdir=paste(output.path),initializationList = ccainit$initlist, priorWeight = 0.9)
    
    ###blurbs start here
    
    e2i_Mn<-matrixToImages(t((myeig2_Mn$eig1)),mask = mask)
    e2i_jac<-matrixToImages(t((myeig2_jac$eig1)),mask = mask)
    e2i_chi<-matrixToImages(t((myeig2_chi$eig1)),mask = mask)
    mask2<-ccainit$initlist[[1]]+ccainit$initlist[[2]]
    
    for (i in 1:2*mynvecs){
      antsImageWrite(ccainit$initlist[[i]],paste(output.path,extension,'ccainit_' ,as.character(i), 'fold', toString(myfold), '_jac.nii.gz',sep=''))
      antsImageWrite((e2i_jac)[[i]],paste(output.path,extension,'sd2eig1_' ,as.character(i), 'fold', toString(myfold), '_jac.nii.gz',sep=''))
      antsImageWrite((e2i_Mn)[[i]],paste(output.path,extension,'sd2eig1_' ,as.character(i), 'fold', toString(myfold), '_Mn.nii.gz',sep=''))
      antsImageWrite((e2i_chi)[[i]],paste(output.path,extension,'sd2eig1_' ,as.character(i), 'fold', toString(myfold), '_chi.nii.gz',sep=''))
      
      }
     
      
    end_time <- Sys.time()
    t1time<-end_time - start_time
    print(t1time)
    
    ###blurbs end here
    
    imgpredtrain_jac<-jac.train %*% (myeig2_jac$eig1)
    imgpredtest_jac<-jac.test %*% (myeig2_jac$eig1)
    imgpredtrain_Mn<-Mn.train %*% (myeig2_Mn$eig1)
    imgpredtest_Mn<-Mn.test %*% (myeig2_Mn$eig1)
    imgpredtrain_chi<-chi.train %*% (myeig2_chi$eig1)
    imgpredtest_chi<-chi.test %*% (myeig2_chi$eig1)
    
    ####start do combo alex
    ncolcombo<-ncol( imgpredtrain_jac)+ncol( imgpredtrain_Mn)+ncol( imgpredtrain_chi)
    
    projs.train <- data.frame(cbind(dist4.train, imgpredtrain_jac , imgpredtrain_Mn ,imgpredtrain_chi)) # column combind the behavior wth the projections
    colnames(projs.train) <- c('Dist_4', paste0('Proj', c(1:ncolcombo)))
    projs.test <- data.frame(cbind(dist4.test, imgpredtest_jac,  imgpredtest_Mn , imgpredtest_chi)) # column combind the behavior wth the projections
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
  ncolcombo<-3*mynvecs
  rows.valid <- c(1:24)
  
  jac.valid <- jac_mat[rows.valid, ]
  Mn.valid <- Mn_mat[rows.valid, ]
  chi.valid <- chi_mat[rows.valid, ]
  
  eig_files_jac <- list.files(path = paste(output.path,sep=''), pattern=paste('*', 'fold', toString(myfold), '_jac.nii.gz', sep=''),full.names = T,recursive = T)
  eig_mat_jac <- imagesToMatrix(eig_files_jac,mask)
  eig_files_Mn <- list.files(path = paste(output.path,sep=''), pattern=paste('*', 'fold', toString(myfold), '_Mn.nii.gz', sep=''),full.names = T,recursive = T)
  eig_mat_Mn <- imagesToMatrix(eig_files_Mn,mask)
  eig_files_chi <- list.files(path = paste(output.path,sep=''), pattern=paste('*', 'fold', toString(myfold), '_chi.nii.gz', sep=''),full.names = T,recursive = T)
  eig_mat_chi <- imagesToMatrix(eig_files_chi,mask)
  
  
  imgmat_Mn_valid <- Mn.valid %*% t(eig_mat_Mn) # [24,numvox] [nvecsx3,numvox]
  imgmat_jac_valid <- jac.valid %*% t(eig_mat_jac)
  imgmat_chi_valid <- chi.valid %*% t(eig_mat_chi)
  
  dist4.valid  <- behavior[rows.valid, 'd4']
  
  projs.valid <- data.frame(cbind(dist4.valid,imgmat_jac_valid,imgmat_Mn_valid,imgmat_chi_valid))
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
  # for (i in 1:numcols){
  #   res1eig<-matrixToImages(jeanat_mang$fusedlist,mask = mask)[[i]] #eig1
  #   antsImageWrite(res1eig,paste(output.path,extension,'JoinEanat' ,as.character(i), '.nii.gz',sep=''))
  # }
  
  mycorsdf_eig2d4<-data.frame(rbind(pcor,corval),row.names=c("pcor","cor"))
  colnames(mycorsdf_eig2d4)<-c('total', paste0('Proj', c(1:ncolcombo)))
  write.csv(mycorsdf_eig2d4, file = paste(output.path ,extension,'fold', toString(myfold), 'd4corsvalidsd2.csv',sep=''))
  
  myperf<-data.frame(rbind(distpred,dist4.valid),row.names=c("d_predicted","d_valid"))
  write.csv(myperf, file = paste(output.path ,extension,'fold', toString(myfold), 'distances4_validsd2.csv',sep=''))
  
} #end for all j labels