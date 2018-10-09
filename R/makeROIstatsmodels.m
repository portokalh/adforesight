mylines1=mybetterroi_modelstats('Hc')
mylines2=mybetterroi_modelstats('PaS')
mylines3=mybetterroi_modelstats('CEnt')
mylines4=mybetterroi_modelstats('fx')
mylines5=mybetterroi_modelstats('M1')
mylines6=mybetterroi_modelstats('M2')

infold=6
allmylines=[mylines1; mylines2; mylines3; mylines4; mylines5; mylines6]

RMSE1=sum((allmylines(:,1:infold)-allmylines(:,25:25+infold-1)).^2,2);
RMSE2=sum((allmylines(:,infold+1:2*infold)-allmylines(:,25+infold:25+2*infold-1)).^2,2);
RMSE3=sum((allmylines(:,2*infold+1:3*infold)-allmylines(:,25+2*infold:25+3*infold-1)).^2,2);
RMSE4=sum((allmylines(:,3*infold+1:4*infold)-allmylines(:,25+3*infold:25+4*infold-1)).^2,2);

mymat=[sqrt(RMSE1/6) sqrt(RMSE2/6) sqrt(RMSE3/6) sqrt(RMSE4/6)]; 
mymat2=[mymat mean(mymat,2) std(mymat')']
dlmwrite('/Users/alex/Documents/GitHub/adforesight/mydata/outdata/ROImodels.csv',mymat2)