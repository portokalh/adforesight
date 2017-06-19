#How good were our original regions as predictors

library( knitr )
library(ANTsR)

# Load in Behavior and Imaging Data
behavior <- read.csv("/Users/omega/Documents/Natalie/SCCAN_tutorial/Eig_zscore/data/All_Behavior.csv")
labled.set <-read.csv("/Users/omega/Documents/Natalie/SCCAN_tutorial/Eig_zscore/data/legendsCHASS2symmetric.csv")
labeled.brain.img <-  antsImageRead("/Users/omega/Documents/Natalie/SCCAN_tutorial/Eig_zscore/data/MDT_labels_chass_symmetric.nii.gz")
mask <-antsImageRead('/Users/omega/Documents/Natalie/SCCAN_tutorial/Eig_zscore/data/MDT_mask_e3.nii')
mang_files<-list.files(path = "/Users/omega/Documents/Natalie/SCCAN_tutorial/Eig_zscore/data", pattern = "T2_to_MDT",full.names = T,recursive = T)
mang_mat <-imagesToMatrix(mang_files,mask)

rows.test <-as.integer(c(2,3,4,9,15,17))
rows.train <-as.integer(c(23,8,18,13,14,19,16,12,11,24,20,7,6,1,22,10,5,21))

mang.train <- mang_mat[rows.train, ]
mang.test <- mang_mat[rows.test, ]
behav.train <- behavior[rows.train, ]
behav.test <-behavior[rows.test, ]

roi.guide <-read.csv('/Users/omega/Documents/MATLAB/connectivity_scripts/ROIguide.csv')
dist4.train <-behav.train[,'d4']
dist4.test  <-behav.test[,'d4']
behav.matrix <-as.matrix(dist4.train)
eanat_region_mang<-sparseDecom2(inmatrix = list(mang.train,behav.matrix), inmask = c(mask, NA), sparseness = c(0.01,1),
                                nvecs = 20, its = 10, cthresh = c(100, 0), mycoption = 0, smooth = 0.01)

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

first <-imageListToMatrix(eanat_region_mang$eig1,mask)
proj.pred <- mang.train %*% t(first)
comb.proj <- data.frame(cbind(dist4.train, proj.pred))
colnames(comb.proj)<-c('Dist_4',paste0( 'Proj', c( 1:ncol(proj.pred ) ))) # insert column names
listpval <-lapply(1:20,function(i){
  y <- paste0('Dist_4')
  x <- paste0('Proj',i)
  LHS <- paste(x)
  form <- as.formula(paste(y, "~",LHS))
  modsum <-summary(lm(formula = form, data = comb.proj))
  pval <- modsum$coefficients[2,4]
})

num.pval <-unlist(listpval)
significantIndices <- which(num.pval < 0.1 )
regions <- reportAnatomy( list(eanat_region_mang$eig1[[significantIndices]] ), mask )
regions
