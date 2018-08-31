
library( knitr )
library(ANTsR)
library(visreg)
library(robustbase)
library(groupdata2)
library(ggplot2)

#reuses fused eigenregions for best performing fold
#Alex Badea 24 august 2018
# Manganese Enhanced MRI   predicts cognitive  performnce
mypath<-'/Volumes/CivmUsers/omega/alex/GitHub/adforesight/'
# Load in Behavior and Imaging Data
mypath <- '/Users/alex/GitHub/adforesight/'
setwd(mypath)


mysp <- 0.05  # 0.01  # 0.05 #0.2 #0.05 # was 0.005 sparseness
mynvecs <- 50 #put back to 50 alex #tested with 10
myell1 <- 1 # make this smaller 0.5, Brian says it is not what i think just switch between l0 and l1
myits<-5 #put back to 5
mysmooth<-0  # was 0
myclus<-250
extension<-paste('JAC', 'sp', toString(mysp), 'vecs', toString(mynvecs), 's', toString(mysmooth),'clus', toString(myclus), sep='') # 'JACsp0p005s0'

output.path <- paste(mypath,'/mydata/outdata2/',extension, '/', sep='') #sd2_projall_noscale/'
#setwd( '/Users/omega/alex/adforesight/' )
source(paste(mypath, '/R/myr2score.R',sep=''))
       
#output.path <- '/Users/omega/alex/adforesight/mydata/outdata/sd2_projall_noscale/'
if (dir.exists(output.path)){ 1} else {dir.create(output.path, recursive=TRUE)}
#if (dir.exists(paste(output.path,'/',extension, '/', sep='')){ 1} else {dir.create(paste(output.path,'/',extension, '/'), recursive=TRUE)}

# Load in Behavior and Imaging Data
behavior <- read.csv('./mydata/All_Behavior.csv')
labled.set <-read.csv('./mydata/legendsCHASS2symmetric.csv')
labeled.brain.img <-  antsImageRead('./mydata/MDT_labels_chass_symmetric.nii.gz')
mask  <- antsImageRead('./mydata/MDT_mask_e3.nii')
mask <- thresholdImage( mask, 0.1, Inf )

# mang_files <- list.files(path = "./mydata/imdata/", pattern = "T2_to_MDT",full.names = T,recursive = T)
# mang_mat <- imagesToMatrix(mang_files,mask)

jac_files <- list.files(path = "./mydata/imdata/", pattern = "jac_to_MDT",full.names = T,recursive = T)
#jac_mat <- imagesToMatrix(jac_files,mask)

#suscept_files <- list.files(path = "./mydata/imdata/", pattern = "X_to_MDT",full.names = T,recursive = T)
#suscept_mat <- imagesToMatrix(suscept_files,mask)

#pick your contrast
mang_files <- jac_files
mang_mat <- imagesToMatrix(mang_files,mask)

rows.test <- as.integer(c(2,3,4,9,15,17))
rows.train <- as.integer(c(23,8,18,13,14,19,16,12,11,24,20,7,6,1,22,10,5,21))

set.seed(1)

#https://cran.r-project.org/web/packages/groupdata2/vignettes/cross-validation_with_groupdata2.html#creating-folds-for-cross-validation  
mygroup <- behavior$genotype[1:24]
myindex <- c(1:24)

mydfb <- data.frame("mysubject_index" = factor(as.integer(myindex)),"mygenotype"=mygroup)
kable(mydfb, align = 'c')

parts <- partition(mydfb, p = c(0.25), id_col = "mysubject_index", cat_col = 'mygenotype')
parts <- partition(mydfb, p = c(0.3), id_col = "mysubject_index", cat_col = 'mygenotype')

test_set <- parts[[1]]
train_set <- parts[[2]]

# %superseeding randomness
test_set$mysubject_index <- rows.test
test_set$mygenotype <- mygroup[rows.test]
train_set$mysubject_index <- rows.train
train_set$mygenotype <- mygroup[rows.train]

# Show test_set
test_set %>% kable()
train_set %>% kable()

train_set <- fold(train_set, k = 4, cat_col = 'mygenotype', id_col = 'mysubject_index')
train_set <- train_set[order(train_set$.folds),]
train_set %>% kable()

# Order by .folds
train_set <- train_set[order(train_set$.folds),]
train_set %>% kable()

set.seed(1)
k<-4



performances <- c()
myBICs  <- c()
myR2score<-c()
myps<-c()


for (fold in 1:k){
 # for (fold in 3){
  gc(verbose = TRUE, reset = FALSE, full = TRUE)
  print('fold:',fold)
  print(fold)

  # Create training set for this iteration
   training_set <- train_set[train_set$.folds != fold, ]
  # training_set %>% kable()
   testing_set <- train_set[train_set$.folds == fold, ]
  # testing_set %>% kable()
  # # Train linear model on training set
  # model <-  lm(model, training_set)
  # # Predict the dependent variable in the testing_set with the trained model
  # predicted <- predict(model, testing_set, allow.new.levels=TRUE)
  # RMSE <- rmse(predicted, testing_set[[dependent]])
  # performances[fold] <- RMSE
  
  
  rows.train <- as.integer(training_set$mysubject_index)
  rows.test <- as.integer(testing_set$mysubject_index)
  mang.train <- mang_mat[rows.train, ]
  mang.test <- mang_mat[rows.test, ]
  
  start_time <- Sys.time()
 
  #alx adds powerit=1, statdir=output.path, robust=1
  #considers scale(mat) cthresh(150) and smooth=0.01
 # eanat_mang<-sparseDecom( inmatrix=mang.train, inmask=mask, nvecs=mynvecs, sparseness=mysp, cthresh=250, its=myits, mycoption=0,verbose=1,statdir=output.path, smooth=mysmooth )
  eanat_mang<-sparseDecom( inmatrix=mang.train, inmask=mask, nvecs=mynvecs, sparseness=mysp, cthresh=myclus, its=myits, mycoption=0,verbose=1,statdir=output.path, smooth=mysmooth )
  
  end_time <- Sys.time()
  t1time<-end_time - start_time
  print(t1time)
  
  behav.train <- behavior[training_set$mysubject_index, ]
  behav.test  <- behavior[testing_set$mysubject_index, ]
  
  gc(verbose = TRUE, reset = FALSE, full = TRUE)
  
  e1l<-list(eanat_mang$eigenanatomyimages)
  
 
  
  
  #alex figure out mask formar
  #graph density range c(1:10)/50
  jeanat_mang <- joinEigenanatomy(mang.train, mask=NA, eanat_mang$eigenanatomyimages, graphdensity=c(0.1,0.2, 0.4,0.6, 0.8), joinMethod='multilevel', verbose=TRUE)
  useeig_mang <- jeanat_mang$fusedlist
  
  #alex prompt
  avgmat_mang <- jeanat_mang$fusedlist
  avgmat_mang <- avgmat_mang/rowSums(abs(avgmat_mang))
  imgmat_train_mang <- (mang.train %*% t(avgmat_mang)  )
  imgmat_test_mang <- (mang.test %*% t(avgmat_mang)  )
  myjvecs<-dim(avgmat_mang)[1]
  
  for (i in 1:myjvecs){
    joint1eig<-matrixToImages(avgmat_mang,mask = mask)[[i]] #eig1
    antsImageWrite(joint1eig,paste(output.path,extension,'JointEigv' ,as.character(i), 'fold', toString(fold), '.nii.gz',sep=''))
  }
  
  
  dist4.train <- behav.train[,'d4']
  dist4.test <- behav.test[,'d4'] 
  
  projs.train <- data.frame(cbind(dist4.train,imgmat_train_mang)) # column combined the behavior wth the projections
  colnames(projs.train) <- c('Dist_4', paste0('Proj', c(1:ncol(imgmat_train_mang)))) # insert column names
  
  projs.test <- data.frame(cbind(dist4.test,imgmat_test_mang)) # column combind the behavior wth the projections
  colnames(projs.test) <- c('Dist_4', paste0('Proj', c(1:ncol(imgmat_test_mang)))) # insert column names
  
  mylm <- lm('Dist_4 ~ .', data=projs.train) # behavior correlation with the number of projections
  
  # step <- stepAIC(mylm, direction="both")
  # step$anova # display results
  
  
  distpred <- predict.lm(mylm, newdata=projs.test) # based on the linear model predict the distances for the same day
  modsum <-summary(lm(distpred~dist4.test))
  r2 <- modsum$adj.r.squared
  my.p <- modsum$coefficients[2,4]
  plot(dist4.test, distpred, xlab = 'Real Swim Distance on Day 4', ylab = 'Predicted Swim Distance on day 4', 
       main='Predicted vs. Real Swim Distance on Day 4', ylim = c(0,1200),xlim = c(0,1200)) # generate plot
  
  
  mymodel<-lm(distpred~dist4.test)
  mytheme <- theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
                   panel.background = element_rect(fill = "white"))
  
  
  RSS <- c(crossprod(mymodel$residuals))
  MSE <- RSS / length(mymodel$residuals)
  RMSE <- sqrt(MSE)
  
  RMSE<-sqrt(mean((distpred - dist4.test)^2)) 
  
  #mymodelset[[fold]]<-mymodel
  performances[fold]<-RMSE2
  myBICs[fold] <- BIC(mymodel)
  myR2score[fold]<-myr2score(distpred,dist4.test)
  myps[fold]<-my.p
  
  myplot<- visreg(mymodel, gg=TRUE) 
  myplot + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
                 panel.background = element_rect(fill = "transparent", colour = NA), 
                 axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                 axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) +
                 ggtitle(paste("RMSE=",formatC(RMSE,digits=2, format="f"), "R2score=",formatC(myR2score[fold],digits=2, format="f"), "  BIC=", formatC(BIC(mymodel),digits=2, format="f")))
  ggsave(paste(output.path,extension,'Mnfold',toString(fold),'.pdf',sep=''), plot = last_plot(), device = 'pdf', 
         scale = 1, width = 4, height = 4, units = c("in"),dpi = 300)
  
  save(mylm, file=paste(output.path , "model", toString(fold), ".Rdata", sep=''))
  #myplot + annotate("text", x = 4, y = 25, label = paste("RMSE=",RMSE))
  
}


###################################
#####    for validation now    ####
###################################

myminfold<-which(performances == min(performances), arr.ind = TRUE)
fold<-myminfold
load(file=paste(output.path , "model", toString(fold), ".Rdata", sep='')) 

#get model for fold 3 or whatever the minimum is
#besteigfiles<-list.files(path=paste(output.path, sep=''), pattern = paste('fold',toString(fold),'.nii.gz',sep=''),full.names = T,recursive = T)
#best_eig_jac_mat <- imagesToMatrix(besteigfiles,mask)
#eigsum<-antsAverageImages(besteigfiles,NORMALIZe=FALSE)
#antsImageWrite(eigsum,paste(output.path,extension,'EigSUM' ,as.character(i), 'fold', toString(fold), '.nii.gz',sep=''))
#myeig<-imagesToMatrix(eigsum,mask)

#read in validation/testing set
rows.valid <- as.integer(test_set$mysubject_index)
mang.valid <- mang_mat[rows.valid, ]

#read joint engen files from best fold
jointeig_files <- list.files(path = paste(output.path,sep=''), pattern=paste('JointEigv','*', 'fold', toString(fold), sep=''),full.names = T,recursive = T)
jointeig_files <- list.files(path = paste(output.path,sep=''), pattern=paste('*', 'fold', toString(fold), '.nii.gz', sep=''),full.names = T,recursive = T)
best_jointeig_mat <- imagesToMatrix(jointeig_files,mask)
imgmat_mang_valid <- mang.valid %*% t(best_jointeig_mat)

#%join eigenanatomy uses this for the best fold
#best_jeanat_mang <- joinEigenanatomy(mang.train, mask=NA, best_eig_jac_mat, graphdensity=0.1, joinMethod='multilevel', verbose=TRUE)
#avgmat_mang <- best_jeanat_mang$fusedlist
#avgmat_mang <- avgmat_mang/rowSums(abs(avgmat_mang))

# i<-1
# imgmat_mang_valid <- (mang.valid %*% t(matrixToImages(eanat_mang$eigenanatomyimages,mask = mask)[[i]])) # load those for fold 2
# e1<-eanat_mang$eigenanatomyimages
# imgmat_mang_valid <- (mang.valid %*% (e1[i,]))

dist4.valid  <- behavior[rows.valid, 'd4']
projs.valid <- data.frame(cbind(dist4.valid,imgmat_mang_valid))
colnames(projs.valid) <- c('Dist_4', paste0('Proj', c(1:ncol(imgmat_mang_valid )))) # insert column names
distpred <- predict.lm(mylm, newdata=projs.valid) 

mymodel<-lm(distpred~dist4.valid)

RSS <- c(crossprod(mymodel$residuals))
MSE <- RSS / length(mymodel$residuals)
RMSE <- sqrt(MSE)

mysummary <-summary(mymodel)
r2pred <- mysummary$adj.r.squared
ppred <- mysummary$coefficients[2,4]


max(behavior$d4[1:24])
#mymodelset[[fold]]<-mymodel
RMSE_valid<-RMSE
BIC_valid <- BIC(mymodel)
R2score_valid<-myr2score(distpred,dist4.valid)

#panel.background = element_rect(fill = "white")
myplot<- visreg(mymodel, gg=TRUE, ylim=c(0,1200),xlim=c(0,1200)) 
myplot + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
               panel.background = element_rect(fill = "transparent", colour = NA), 
               axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
               axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) +
    #xlim(0,1200)+ylim(0,1200)+coord_cartesian(xlim = c(1200, 1200),ylim = c(1200, 1200)) + coord_equal()+
  ggtitle(paste("RMSE=",formatC(RMSE_valid,digits=2, format="f"), "R2score=",formatC(R2score_valid,digits=2, format="f"), "  BIC=", formatC(BIC_valid,digits=2, format="f"), 
                " R2=", formatC(r2pred,digits=2, format="f"), " p= ", formatC(ppred,digits=4, format="f")))

ggsave(paste(output.path,extension,'MnValidationSet',toString(fold),'.pdf',sep=''), plot = last_plot(), device = 'pdf', 
       scale = 1, width = 4, height = 4, units = c("in"),dpi = 300)

numcols<-dim(imgmat_mang_valid)[2]
rd4 <- t(t(dist4.valid)[rep(1,c(numcols)),])
pcor<-c(numcols)
corval<-c(numcols)

cor(rd4,imgmat_mang_valid)
for (i in 1:numcols) {
  mypcor<-cor.test(t(dist4.valid),t(imgmat_mang_valid[,i]))
  pcor[i]<-mypcor$p.value
  corval[i]<-mypcor$estimate
}

# for (i in 1:numcols){
#   res1eig<-matrixToImages(jeanat_mang$fusedlist,mask = mask)[[i]] #eig1
#   antsImageWrite(res1eig,paste(output.path,extension,'JoinEanat' ,as.character(i), '.nii.gz',sep=''))
# }

mycorsdf_eig2d4<-data.frame(rbind(pcor,corval),row.names=c("pval","cor"))
write.csv(mycorsdf_eig2d4, file = paste(output.path ,extension,'eig2d4cors.csv',sep=''))
#redo fold min or save models

# reportAnatomy <- function( eigenImageList, mask, weight = 0.3 )
#   {
#   sccanAalLabels <- c()
#   for( eigenImage in eigenImageList )
#   {
#     nonZeroIndices<- abs( eigenImage[mask == 1] ) > 0
#     sccanAalLabels <- append( sccanAalLabels, aalImage[mask == 1][nonZeroIndices] )
#   }
# }
#   reportedAnatomy <- reportAnatomy( SccanImages, mask )