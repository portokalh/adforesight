#Correlation Between Imaging and Swim Distance

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

label.ind <- c(19,20,42,47,51,54,57,59,62,63,64,65,73,81,83,91,119,121,123)

#Both (Left and Right) ROIS


for ( j in 1:length(label.ind)) {
  labels<-c(label.ind[j],label.ind[j]+1000)
  initmat<-matrix( rep(0,sum(mask==1)*length(labels)), nrow=length(labels) )
  output.file <- paste(output.path, as.character(labled.set$Abbreviation[label.ind[j]]),'LR_all.png',sep = "_")
  #fill the matrix with the aal region locations
  for ( i in 1:length(labels) ) {
    vec<-( labeled.brain.img[ mask == 1 ] == labels[i] )
    initmat[i,]<-as.numeric( vec )
  }
  newmat <- as.matrix(colSums(initmat))
  newmat <- t(newmat)
  ccainit<-initializeEigenanatomy(newmat,mask)
  
  dist4.train <-behav.train[,'d4']
  dist4.test  <-behav.test[,'d4']
  behav.matrix <-as.matrix(dist4.train)
  eanat_region_mang<-sparseDecom2(inmatrix = list(mang.train,behav.matrix), inmask = c(mask, NA), sparseness = c(0.01,1),
                                  nvecs = 1, its = 10, cthresh = c(100, 0), mycoption = 0, smooth = 0.01,
                                  initializationList = ccainit$initlist,priorWeight = 0.5)
  
  imgmat_mang<-imageListToMatrix( eanat_region_mang$eig1 , mask )
  imgpredtrain_mang<-mang.train %*% t(imgmat_mang)
  imgpredtest_mang<-mang.test %*% t(imgmat_mang)
 
  projs.train <- data.frame(cbind(dist4.train, imgpredtrain_mang)) # column combind the behavior wth the projections
  colnames(projs.train) <- c('Dist_4','Mang') # insert column names
  projs.test <- data.frame(cbind(dist4.test, imgpredtest_mang)) # column combind the behavior wth the projections
  colnames(projs.test) <- c('Dist_4','Mang') # insert column names
  mylm <- lm('Dist_4 ~ .', data=projs.train) # behavior correlation with the number of projections
  distpred <- predict.lm(mylm, newdata=projs.test) # based on the linear model predict the distances for the same day
  modsum <-summary(lm(distpred~dist4.test))
  r2 <- modsum$adj.r.squared
  my.p <- modsum$coefficients[2,4]
  png(filename = output.file, width = 556, height = 579,units = "px")
  plot(imaging,cognition,ylab = 'Projection on Swim Distance', xlab = 'All Modalties ',
       main='Correlation Between Cognition vs Imaging') # generate plot
  abline(lm(myform))
  rp = vector('expression',2)
  rp[1] = substitute(expression(italic(R)^2 == MYVALUE), list(MYVALUE = format(r2, digits = 3)))[2]
  rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
  legend('top', legend = rp, bty = 'n')
  dev.off()
}
