%% Script to process wilting/critical point from 20yr NNsm data
 % Author: Qing He
 
 data_prefix = '/home/xyw/matlabcode/nnsm_layer1/';
 
 yrlist = 2002:2023;
 thetaf_allyr = nan(406,964,length(yrlist));
 thetac_allyr = nan(406,964,length(yrlist));

 
 for yr = 1:length(yrlist)
     disp(yrlist(yr))
     load([data_prefix,'NNsm_mdata_tauL_f1_',...
           num2str(yrlist(yr)),'.mat'],'thetaf_allm','psm_allm');
       
     thetaf_allyr(:,:,yr) =   thetaf_allm;
     thetac_allyr(:,:,yr) =   psm_allm;
     
 end
 
 
 thetaf_allyr_median = nanmedian(thetaf_allyr,3);
 thetac_allyr_median = nanmedian(thetac_allyr,3);
 
 
 save('NNSM_22year_soil_thresholds_f1.mat', 'thetaf_allyr', 'thetac_allyr','thetaf_allyr_median','thetac_allyr_median','-v7.3')

                        
 %% process based on raw estimation
 lat = nan(406,964);
 data_prefix = '/home/xyw/matlabcode/nnsm_layer1/soil_memory/'; 
 yrlist = 2002:2023;
 
 for yr = 1:length(yrlist)
     disp(yrlist(yr))
     load([data_prefix,'NNsm_tauL_f1_',...
           num2str(yrlist(yr)),'.mat'],'lamuda','theta0','thetaf','psm','R2');
 
     r = [];
     for i = 1:size(theta0,1)
         if isnan(theta0(i,1))
             break
         end
         %if theta0(i,3) < 0.05*SMrange(theta0(i,1),theta0(i,2)) || theta0(i,3) < 2*0.04 % SMAP has an error (ubRMSD) <= 0.04
         if theta0(i,3) < 2*0.04 %3*0.04 %0.04 % SMAP has an error (ubRMSD) <= 0.04
             r = [r; i];
         end
     end
     lamuda(r,:) = [];
     thetaf(r,:) = [];
     theta0(r,:) = [];
     psm(r,:) = [];
     R2(r,:) = [];
     
     thetaf_all = nan(size(lat,1),size(lat,2),100);
     psm_all = nan(size(lat,1),size(lat,2),100);
     theta0_all = nan(size(lat,1),size(lat,2),100);
     R2_all = nan(size(lat,1),size(lat,2),100);
     count_all = nan(size(lat,1),size(lat,2));
     
     for i = 1:size(lamuda,1)
         %disp(i)
        if isnan(lamuda(i,1))
            break
        end
        if i == 1
            count = 1;
        elseif lamuda(i,1) ~= lamuda(i-1,1) || lamuda(i,2) ~= lamuda(i-1,2)
            count = 1;
        else
            count = count + 1;
        end
        
        thetaf_all(lamuda(i,1),lamuda(i,2),count) = thetaf(i,3);
        psm_all(lamuda(i,1),lamuda(i,2),count) = psm(i,3);
        theta0_all(lamuda(i,1),lamuda(i,2),count) = theta0(i,3);
        R2_all(lamuda(i,1),lamuda(i,2),count) = R2(i,3);
        count_all(lamuda(i,1),lamuda(i,2)) = count;
     end
     
     thetaf_all(R2_all<=0.7) = NaN;
     theta0_all(R2_all<=0.7) = NaN;
     psm_all(R2_all<=0.7) = NaN;
     
     thetaf_all(count_all<=10) = NaN;
     theta0_all(count_all<=10) = NaN;
     psm_all(count_all<=10) = NaN;
     
     fname = ['NNsm_tauL_fullpro_f3_',num2str(yrlist(yr)),'.mat'];
     save(fname, 'thetaf_all','psm_all','R2_all','count_all','theta0_all')
     
 end
     
 %%
 %data_prefix = '/Volumes/heqing_drive/projects_backup/nmp_sm/NNsm_analyses/';
 data_prefix = '/home/xyw/matlabcode/';
 yrlist = 2002:2023;
  
 thetaf_allyr = [];
 psm_allyr = [];
 theta0_allyr = [];
 for yr = 1:length(yrlist)
     disp(yrlist(yr))
     load([data_prefix,'NNsm_tauL_fullpro_f3_',...
           num2str(yrlist(yr)),'.mat'],'thetaf_all','psm_all','theta0_all');
             
     thetaf_allyr = cat(3,thetaf_allyr,thetaf_all);
     psm_allyr = cat(3,psm_allyr,psm_all);
     theta0_allyr = cat(3,theta0_allyr,theta0_all);
 end
 
 thetaf_std = std(thetaf_allyr,0,3,'omitnan');
 thetaf_cv = thetaf_std./nanmean(thetaf_allyr,3);
 thetaf_5prc = prctile(thetaf_allyr, 5, 3);
 
 psm_std = std(psm_allyr,0,3,'omitnan');
 psm_cv = psm_std./nanmean(psm_allyr,3);
 psm_95prc = prctile(psm_allyr, 95, 3);
 
 theta0_std = std(theta0_allyr,0,3,'omitnan');
 theta0_cv = theta0_std./nanmean(theta0_allyr,3);
 theta0_95prc = prctile(theta0_allyr, 95, 3);

 
 save NNsm_thetafc_statistics_f3.mat thetaf_std thetaf_cv thetaf_5prc ...
                                  psm_std psm_cv psm_95prc ...
                                  theta0_std theta0_cv theta0_95prc
 %% replace thetaf_median with thetaf_5prc over humid regions, and 
  % thetac_median with thetac_95prc over arid regions
  % AI = Rn/(lamuda*P); AI <0.2 for arid regions; AI>=0.65 for humid region

  load('./NNSM_22year_soil_thresholds_f3.mat','thetaf_allyr_median','thetac_allyr_median')
  
  thetaf = thetaf_allyr_median;
  thetac = thetac_allyr_median;

load('AI_mean_era5.mat')
thetaf(AI>=0.65)=thetaf_5prc(AI>=0.65);
thetac(AI<0.2)=psm_95prc(AI<0.2);
 save ./NNsm_soil_thresholds_f3.mat thetaf thetac
 
 
 