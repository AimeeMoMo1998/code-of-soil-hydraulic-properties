%% Script to process wilting/critical point from 20yr NNsm data
 % author Yawei Xu

load('vwc_mask.mat')
vwc_mask(vwc_mask>5)=nan;
vwc_mask(vwc_mask<=5)=1;

%% compact all file
thetaf_allm_allyear_f1=nan(406*964,22);thetaf_allm_allyear_f3=nan(406*964,22);
psm_allm_allyear_f1=nan(406*964,22);psm_allm_allyear_f3=nan(406*964,22);

for year=2002:1:2023
    load(strcat('D:\科研项目\drydown\数据\NNsm_mdata_tauL_f1_',string(year),'.mat'));   
    thetaf_allm_allyear_f1(:,year-2001)=reshape(thetaf_allm,406*964,1);%boxplot(thetaf_allm),title(string(year)),pause(1);
    psm_allm_allyear_f1(:,year-2001)=reshape(psm_allm,406*964,1);%boxplot(psm_allm),title(string(year)),pause(1);   
   
end

%% remove outliers
X=thetaf_allm_allyear_f1;
Y = prctile(X,[25 50 75],1);
IQR=Y(3,:)-Y(1,:);
Upperlimit=Y(3,:)+1.5.*IQR;
Lowerlimit=Y(1,:)-1.5.*IQR;
thetaf_rmoutlier=X;
thetaf_rmoutlier(thetaf_rmoutlier>Upperlimit)=nan;
thetaf_rmoutlier(thetaf_rmoutlier<Lowerlimit)=nan;
r_thetaf_rmoutlier_f1=reshape(thetaf_rmoutlier,406,964,22);

X=psm_allm_allyear_f1;
Y = prctile(X,[25 50 75],1);
IQR=Y(3,:)-Y(1,:);
Upperlimit=Y(3,:)+1.5.*IQR;
Lowerlimit=Y(1,:)-1.5.*IQR;
psm_rmoutlier=psm_allm_allyear_f1;
psm_rmoutlier(psm_rmoutlier>Upperlimit)=nan;
psm_rmoutlier(psm_rmoutlier<Lowerlimit)=nan;
r_psm_rmoutlier_f1=reshape(psm_rmoutlier,406,964,22);

%% global median map
load('mask1.mat','mask');

% global median
thetaf_median_f1=reshape(median(r_thetaf_rmoutlier_f1,3,'omitnan'),406,964);% 22year
psm_median_f1=reshape(median(r_psm_rmoutlier_f1,3,'omitnan'),406,964);% 22year

nnsm_pwp=thetaf_median_f1.*mask.*vwc_mask;
nnsm_cp=psm_median_f1.*mask.*vwc_mask;

load('AI_mean_era5.mat') 
AI_psm_mask=find((AI.*0.1)<0.2);
AI_thetaf_mask=find((AI.*0.1)>=0.65);

load('AI_correction.mat','CP_AI','PWP_AI'); 

replace_data=CP_AI;
nnsm_cp(AI_psm_mask)=replace_data(AI_psm_mask);

i=10;
replace_data=PWP_AI;
nnsm_pwp(AI_thetaf_mask)=replace_data(AI_thetaf_mask);

save('multiyear_nnsm_cp_pwp.mat','nnsm_cp','nnsm_pwp');
