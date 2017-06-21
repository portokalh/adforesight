#Predicting Cognitive Battery

library( knitr )
library(ANTsR)

# Load in Behavior and Imaging Data
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

roi.guide <-read.csv('/Users/omega/alex/adforesight/mydata/ROIguide.csv')
reportAnatomy<-function( eigIn, maskIn, wt=0.3)
{
  ccaanat<-list()
  for ( img in eigIn ) {
    nzind<-abs(img[ maskIn == 1 ]) > 0
    aalvals<-labeled.brain.img[ maskIn == 1 ][ nzind ]
    ccaanat<-lappend( ccaanat, aalvals )
  }
  ccaanat<-unlist( ccaanat )
  anatcount<-hist(ccaanat,breaks=0:1167, plot = F)$count
  anatcount[ anatcount < wt*max(anatcount) ]<-0
  anatcount<-which( anatcount > 0 )
  rois <- match(anatcount,roi.guide$Label_Num)
  return( toString(roi.guide$Abbreviation[rois] ) )
}

# Try to predict all the demographic variability from the imaging data. We use
# `mycoption 0` to try to reduce correlation in low-dimensional space. This
# enforces a new SCCAN constraint (not previously reported). (Nick: We've been
# using mycoption = 0 this whole time.)
behav.train.mat<- as.matrix( behav.train[,c(2:16)] )
behav.test.mat<- as.matrix( behav.test[,c(2:16)] )
cognitiveSccanResult <- sparseDecom2(inmatrix = list (mang.train,behav.train.mat),its = 10, mycoption = 0, sparseness = c( 0.03, 0.05 ), nvecs = 11,
                                     inmask = c( mask, NA ), cthresh = c( 100, 0 ), smooth = 0.01 )
imgmat_mang<-imageListToMatrix( cognitiveSccanResult$eig1 , mask )
for ( bestpred in 1:nrow(imgmat_mang)) {
  imgpredtrain_mang<-mang.train %*% t(imgmat_mang)
  imgpredtest_mang<-mang.test %*% t(imgmat_mang)
  cogpredtrain_mang <- behav.train.mat %*% as.matrix(cognitiveSccanResult$eig2)[,bestpred]
  cogpredtest_mang <- behav.test.mat %*% as.matrix(cognitiveSccanResult$eig2)[,bestpred]
  gvars<-paste("Proj",c(1:nrow(imgmat_mang)),sep='',collapse='+')
  projs.train <- data.frame(cbind(cogpredtrain_mang,imgpredtrain_mang)) # column combind the behavior wth the projections
  colnames(projs.train) <- c('cognitive',paste0( 'Proj', c( 1:ncol(imgpredtrain_mang ) ))) # insert column names
  projs.test <- data.frame(cbind(cogpredtest_mang, imgpredtest_mang)) # column combind the behavior wth the projections
  colnames(projs.test) <- c('cognitive',paste0( 'Proj', c( 1:ncol(imgpredtest_mang ) ))) # insert column names
  myform<-as.formula( paste("cognitive~",gvars,sep='') )
  mylm <- lm(myform, data=projs.train) # behavior correlation with the number of projections
  distpred <- bigLMStats( mylm)
  cat(paste("Eig",bestpred,"is related to:\n"))
  mycog<-colnames(behav.train.mat)[ abs(cognitiveSccanResult$eig2[,bestpred]) > 0 ]
  cat( mycog )
  cat("\nwith weights\n")
  cat( abs(cognitiveSccanResult$eig2[,bestpred])[ abs(cognitiveSccanResult$eig2[,bestpred]) > 0 ])
  cat(paste("\nwith predictive correlation:",
            cor( cogpredtest_mang,predict(mylm,newdata=projs.test))))
  cat("\nAnatomy:")
  for ( x in which.min(p.adjust(distpred$beta.pval)) )  {
    myanat<-reportAnatomy( list( cognitiveSccanResult$eig1[[x]]) , mask , 0.3 )
    cat(myanat)
    # if ( render ) {
    #   vizimg<-antsImageClone( batt$eig1[[x]] )
    #   ImageMath(3,vizimg,'abs',vizimg)
    #   brain<-renderSurfaceFunction( surfimg =list( aalimg ) , alphasurf=0.1 ,
    #                                 funcimg=list(vizimg), smoothsval = 1.5 )
    #   id<-par3d("userMatrix")
    #   rid<-rotate3d(  id , -pi/2, 1, 0, 0 )
    #   rid2<-rotate3d( id ,  pi/2, 0, 0, 1 )
    #   rid3<-rotate3d( id , -pi/2, 0, 0, 1 )
    #   par3d(userMatrix = id )
    #   ofn<-paste(rootdir,'/figures/battery',bestpred,sep='')
    #   nfn[ bestpred ]<-paste(ofn,'.png',sep='')
    #   cognames[ bestpred ]<-paste(mycog,collapse='+')
    #   dd<-make3ViewPNG(  rid, id, rid2, ofn )
    #   par3d(userMatrix = id )
    # }
    cat("\n")
  }
  cat("\n")
}
