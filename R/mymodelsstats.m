%Hc
mypath='/Users/alex/Documents/GitHub/adforesight/mydata/outdata/VBA_jac_posneg_pred/'
part1=csvread([mypath, '/VBA_jac_posneg_preddistances4_pv_fold1.csv'],1,1)
part2=csvread([mypath, '/VBA_jac_posneg_preddistances4_pv_fold2.csv'],1,1)
part3=csvread([mypath, '/VBA_jac_posneg_preddistances4_pv_fold3.csv'],1,1)
part4=csvread([mypath, '/VBA_jac_posneg_preddistances4_pv_fold4.csv'],1,1)
Hcline=[part1(1,:) part2(1,:) part3(1,:) part4(1,:) part1(2,:) part2(2,:) part3(2,:) part4(2,:)]

mypath='/Users/alex/Documents/GitHub/adforesight/mydata/outdata/VBA_jac_posneg_pred/';
part1=csvread([mypath,'VBA_jac_posneg_preddistances4_pv_fold1.csv'],1,1)
part2=csvread([mypath,'VBA_jac_posneg_preddistances4_pv_fold2.csv'],1,1)
part3=csvread([mypath,'VBA_jac_posneg_preddistances4_pv_fold3.csv'],1,1)
part4=csvread([mypath,'VBA_jac_posneg_preddistances4_pv_fold4.csv'],1,1)
Hc_vol=[part1(1,:) part2(1,:) part3(1,:) part4(1,:) part1(2,:) part2(2,:) part3(2,:) part4(2,:)]

mypath='/Users/alex/Documents/GitHub/adforesight/mydata/outdata/VBA_Mn_posneg_pred/'
part1=csvread([mypath, '/VBA_Mn_posneg_preddistances4_pv_fold1.csv'],1 ,1);



mypath='/Users/alex/Documents/GitHub/adforesight/mydata/outdata/combo_ROI_predsp0.05smooth0.1/Hc/'
part1=csvread([mypath,'/combo_ROI_preddistances4_pv_fold1.csv'],1,1) ;
part2=csvread([mypath,'/combo_ROI_preddistances4_pv_fold2.csv'],1,1);
part3=csvread([mypath,'/combo_ROI_preddistances4_pv_fold3.csv'],1,1);
part4=csvread([mypath,'/combo_ROI_preddistances4_pv_fold4.csv'],1,1)
Hc_combo=[part1(1,:) part2(1,:) part3(1,:) part4(1,:) part1(2,:) part2(2,:) part3(2,:) part4(2,:)]
