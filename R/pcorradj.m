myprefix={'Hc', 'PaS', 'CEnt', 'fx', 'M1', 'M2'}
mycontrast={'jac', 'Mn' , 'chi'}


for i=1:numel(myprefix)
    %for i=1:2
    i
    prefix=char(myprefix{i})
    for j=1:3
        contrast=mycontrast{j};
    
        mypath=['/Users/alex/Documents/GitHub/adforesight/mydata/outdata/single_ROI_predsp0.05smooth0.1/', contrast, '/', prefix, '/single_ROI_predfold4d4corsvalidsd2.csv'];
        res=csvread(mypath,1,1);
        myp((i-1)*4+j) = res(1,1); mycorr((i-1)*4+j)= res(2,1);
    end
    for j=4
        mypath=['/Users/alex/Documents/GitHub/adforesight/mydata/outdata/combo_ROI_predsp0.05smooth0.1/' prefix '/combo_ROI_predfold4d4corsvalidsd2.csv'];  
        res=csvread(mypath,1,1);
        myp((i-1)*4+j) = res(1,1); mycorr((i-1)*4+j)= res(2,1);
    end
end

r2=mycorr.^2;
adjR2=1-(1-r2)*23/(24-2);
mymat=[myp ;mycorr ;r2 ;adjR2]';
dlmwrite('/Users/alex/Documents/GitHub/adforesight/mydata/outdata/ROImodelsWHOLE.csv',[myp ;mycorr ;r2 ;adjR2]' )