%% Script to estimate tauL of reanalysis data
 % Author: Qinghe, Mar 2022
 % Datasets are: NNsm
 %f=3;
 
 %% post-processing of the multidata results
 load('/mnt/WCL-S920B_home/qinghe/reanaly_compar/latlonP.mat');

 yrlist={'2002','2003','2004','2005','2006','2007','2008','2009'...
     '2010','2011','2012','2013','2014','2015','2016','2017','2018','2019','2020','2021'};
 
 dpsoil1 = 30;  %mm
 tot_yrs = 1;
  
 lamuda_allm = cell(1,1);
 thetaf_allm = cell(1,1);
 psm_allm = cell(1,1);
 R2_allm  = cell(1,1);
 dd_st_allm = cell(1,1);
 dd_end_allm = cell(1,1);
 tau_allm = cell(1,1);
 count_allm = cell(1,1);
 tot_ET_II_allm = cell(1,1);
 
 for d = 1:length(yrlist)
     load(['./NNsm_tauL_f1_',yrlist{d},'.mat']);
     %remove fitting when R2<0.7
     ind = find(R2(:,3)<0.7); %0.9);
     lamuda(ind,:) = [];
     psm(ind,:) = [];
     R2(ind,:) = [];
     theta0(ind,:) = [];
     dd_dates(ind,:) = [];
     thetaf(ind,:) = [];

     %remove fitting caused by noise
     r = [];
     for i = 1:size(theta0,1)
         disp([yrlist{d},', loop1: ',num2str(i)])
        if isnan(theta0(i,1))
            break
        end
        %if theta0(i,3) < 0.05*SMrange(theta0(i,1),theta0(i,2)) || theta0(i,3) < 2*0.04 % SMAP has an error (ubRMSD) <= 0.04
        if theta0(i,3) < 2*0.04 %3*0.04 %0.04 % SMAP has an error (ubRMSD) <= 0.04
            r = [r; i];
        end
     end
     
     theta0(r,:) = [];
     lamuda(r,:) = [];
     thetaf(r,:) = [];
     psm(r,:) = [];
     R2(r,:) = [];
     dd_dates(r,:) = [];
    
     %reshape 
     lamuda_all = nan(size(lat,1),size(lat,2),100);
     theta0_all = nan(size(lat,1),size(lat,2),100);
     thetaf_all = nan(size(lat,1),size(lat,2),100);
     psm_all = nan(size(lat,1),size(lat,2),100);
     R2_all = nan(size(lat,1),size(lat,2),100);
     dd_st_all = nan(size(lat,1),size(lat,2),100);
     dd_end_all = nan(size(lat,1),size(lat,2),100);
     count_all = nan(size(lat,1),size(lat,2));

     for i = 1:size(lamuda,1)
         disp([yrlist{d},', loop2: ',num2str(i)])
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
    
         lamuda_all(lamuda(i,1),lamuda(i,2),count) = lamuda(i,3);
         theta0_all(lamuda(i,1),lamuda(i,2),count) = theta0(i,3);
         thetaf_all(lamuda(i,1),lamuda(i,2),count) = thetaf(i,3);
         psm_all(lamuda(i,1),lamuda(i,2),count) = psm(i,3);
         R2_all(lamuda(i,1),lamuda(i,2),count) = R2(i,3);
         dd_st_all(lamuda(i,1),lamuda(i,2),count) = dd_dates(i,3);
         dd_end_all(lamuda(i,1),lamuda(i,2),count) = dd_dates(i,4);
         count_all(lamuda(i,1),lamuda(i,2)) = count;
     end

     tau_all = 1./lamuda_all;

     %remove fitting when drydown events <3
     minDrydowns = 3; %4;
     lamuda_all(repmat(count_all,[1 1 100])<minDrydowns) = NaN; % filter out regions with < minDrydowns drydown events
     theta0_all(repmat(count_all,[1 1 100])<minDrydowns) = NaN; 
     thetaf_all(repmat(count_all,[1 1 100])<minDrydowns) = NaN; 
     psm_all(repmat(count_all,[1 1 100])<minDrydowns) = NaN; 
     R2_all(repmat(count_all,[1 1 100])<minDrydowns) = NaN;
     dd_st_all(repmat(count_all,[1 1 100])<minDrydowns) = NaN; 
     dd_end_all(repmat(count_all,[1 1 100])<minDrydowns) = NaN; 
     tau_all(repmat(count_all,[1 1 100])<minDrydowns) = NaN; 
     
     lamuda_allm = nanmedian(lamuda_all,3);
     thetaf_allm = nanmedian(thetaf_all,3);
     psm_allm = nanmedian(psm_all,3);
     R2_allm = nanmedian(R2_all,3);
     dd_st_allm = nanmedian(dd_st_all,3);
     dd_end_allm = nanmedian(dd_end_all,3);
     tau_allm = nanmedian(tau_all,3);
     count_allm = count_all;
     
     ET_II_all = dpsoil1.*theta0_all.*...
                 (1-exp(-(dd_end_all-dd_st_all+1)./tau_all));
     tot_ET_II_all = nansum(ET_II_all,3)/tot_yrs;   %mm/yr
     tot_ET_II_all(tot_ET_II_all==0) = NaN;
     tot_ET_II_allm = tot_ET_II_all;
     
%      filename = ['/home/xyw/matlabcode/nnsm_layer1/NNsm_mdata_tauL_f1_',yrlist{d},'.mat'];
%      save(filename,'lamuda_allm', 'thetaf_allm', ...
%                    'psm_allm','R2_allm',...
%                    'dd_st_allm','dd_end_allm','tau_allm',...
%                    'count_allm','tot_ET_II_allm','-v7.3');
      filename = ['/home/xyw/matlabcode/nnsm_layer1/NNsm_tot_ET_II_all_f1_',yrlist{d},'.mat'];
     save(filename,'tot_ET_II_allm','-v7.3');          
 end  