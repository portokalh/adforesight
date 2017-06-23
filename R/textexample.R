data(aal,package='ANTsR')
rootdir='/Users/alex/GitHub/sccantutorial/sccanTutorial-gh-pages/data/'
gfnl<-list.files(path=rootdir, pattern = glob2rx("pbac*mha"),
                 full.names = T,recursive = T)
ptrainimg<-as.matrix(antsImageRead(gfnl[2],2))
ptestimg<-as.matrix(antsImageRead(gfnl[1],2))
gfnl<-list.files(path=rootdir, pattern = "gmask.nii.gz",
                 full.names = T,recursive = T)
mask<-antsImageRead( gfnl[1], 3 )
afnl<-list.files(path=rootdir, pattern = "aal.nii.gz",
                 full.names = T,recursive = T)
aalimg<-antsImageRead( afnl[1], 3 )
f1<-list.files(path =rootdir, pattern = "pbac_train_cog.csv",
               recursive=TRUE, full.names = TRUE, include.dirs=TRUE )
f2<-list.files(path = rootdir, pattern = "pbac_test_cog.csv",
               recursive=TRUE, full.names = TRUE )
ptraincog<-read.csv(f1)
ptestcog<-read.csv(f2)

agemat<-matrix( ptraincog$age, ncol=1)
paramsearch<-c(1:10)/(-100.0)
paramsearchcorrs<-rep(0,length(paramsearch))
paramsearchpreds<-rep(0,length(paramsearch))
ct<-1
for ( sp in paramsearch ) {
  ageresult<-sparseDecom2( inmatrix=list(ptrainimg,agemat), its=8, mycoption=1,
                           sparseness=c(sp,0.9), inmask=c(mask,NA),nvecs=2, cthresh=c(50,0))
  # convert output images to matrix so we can validate in test data
  ccamat<-imageListToMatrix( ageresult$eig1, mask )
  agepred<-ptrainimg %*% t(ccamat)
  paramsearchcorrs[ct]<-cor( agepred[,1],  ptraincog$age )
  agepred<-ptestimg %*% t(ccamat)
  paramsearchpreds[ct]<-cor( agepred[,1],  ptestcog$age )
  ct<-ct+1
}
mydf<-data.frame( sparseness=paramsearch, trainCorrs=paramsearchcorrs,
                  testCorrs=paramsearchpreds )
mdl1<-lm( trainCorrs ~ stats::poly(sparseness,4), data=mydf )
mdl2<-lm( testCorrs ~ stats::poly(sparseness,4) , data=mydf )
visreg(mdl1)