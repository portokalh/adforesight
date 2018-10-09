function mylines=mybetterroi_modelstats(prefix)
%prefix='Hc'
res1_jac=[0, 0]'
res1_T1=[0, 0]'
res1_qsm=[0, 0]'
res1_combo=[0, 0]'

for myfold=1:4
  myfold  
path_jac=['/Users/alex/Documents/GitHub/adforesight/mydata/outdata/single_ROI_predsp0.05smooth0.1/jac/', prefix ,'/single_ROI_preddistances4_pv_fold' num2str(myfold) '.csv'];
path_T1  = ['/Users/alex/Documents/GitHub/adforesight/mydata/outdata/single_ROI_predsp0.05smooth0.1/Mn/', prefix, '/single_ROI_preddistances4_pv_fold' num2str(myfold) '.csv']; 
path_qsm=['/Users/alex/Documents/GitHub/adforesight/mydata/outdata/single_ROI_predsp0.05smooth0.1/chi/', prefix,'/single_ROI_preddistances4_pv_fold' num2str(myfold) '.csv'];
path_combo=['/Users/alex/Documents/GitHub/adforesight/mydata/outdata/combo_ROI_predsp0.05smooth0.1/' prefix '/combo_ROI_preddistances4_pv_fold' num2str(myfold) '.csv']; 

res1_jac=[res1_jac  csvread(path_jac,1,1)]
res1_T1=[res1_T1 csvread(path_T1,1,1)];
res1_qsm=[res1_qsm csvread(path_qsm,1,1)];
res1_combo=[res1_combo csvread(path_combo,1,1)];

end

res1_jac=res1_jac(:,2:end)
res1_T1=res1_T1(:,2:end)
res1_qsm=res1_qsm(:,2:end)
res1_combo=res1_combo(:,2:end)

mylines=[res1_jac(1,:) res1_jac(2,:); res1_T1(1,:) res1_T1(2,:);...
    res1_qsm(1,:) res1_qsm(2,:);res1_combo(1,:) res1_combo(2,:)];
