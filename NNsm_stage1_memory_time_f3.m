cd('/home/xyw/matlabcode/nnsm_layer1/');
clc
clear
load('/mnt/WCL-S920B_home/qinghe/reanaly_compar/latlonP.mat');
% yrlist=...
%     {'2002','2003','2004','2005','2006',...
%     '2007','2008' ,'2009','2010','2011'...
%      '2012','2013','2014','2015','2016'...
%      '2017','2018','2019'};
yrlist={'2020'};
monthlist={'01','02','03','04','05','06','07','08','09','10','11','12'};
% GPM_prec=nan(406,964,365*18+4);


dayindex=0;
for yr=1:length(yrlist)
%     %% load GPM data
%     if mod(str2num(yrlist{yr}),4)==0
%         daymax=366;
%     else
%         daymax=365;
%     end
%     
%     for month=1:1:12
%         GPM_prec_files=dir(strcat('/home/xyw/GPM-IMERG-HDF/GPM to EASE/',yrlist{yr},monthlist{month},'/*.mat'));
%         for day=1:1:size(GPM_prec_files,1)
%             dayindex=dayindex+1;
%             load(strcat(GPM_prec_files(day).folder,'/',GPM_prec_files(day).name));
%             GPM_prec(:,:,dayindex)=rain;
%         end
%     end
    load(strcat('/home/xyw/matlabcode/GPM_process_AFF/','GPM_EASE_Daily_6AM_Local_36km_',yrlist{yr},'01_',yrlist{yr},'12.mat'));
    GPM_prec=GPM_EASE_Daily_6AM_Local;
    %% Stage-I memory time
     load(['/home/xyw/matlabcode/nnsm_layer1/NNsm_year',yrlist{yr},'.mat']);
     SSM_Data = nnsm_year;
     SSM_Data =permute(SSM_Data,[2 1 3]);

      SSM_Data_us = single(nan(size(SSM_Data)));
     f = 1/3;
     for ii = 2:1/f:(size(SSM_Data,3)-1)
        %disp(ii)
        sm_running = SSM_Data(:,:,ii);
        temp = SSM_Data(:,:,ii+1);
        sm_running(~isnan(temp)) = temp(~isnan(temp));
        temp = SSM_Data(:,:,ii-1);
        sm_running(~isnan(temp)) = temp(~isnan(temp));
        SSM_Data_us(:,:,ii) = sm_running;
     end
     
     sm_sum = nan(size(SSM_Data_us,1),size(SSM_Data_us,2));
     for i = 1:size(lat,1)
         for j = 1:size(lat,2)
              sm = reshape(SSM_Data_us(i,j,:),1,size(SSM_Data_us,3));
              sm(sm<0.04) = nan;
              sm(sm>0.7) = nan;

              sm(isnan(sm)) = [];

              if ~isempty(sm)
                 psm = sm(2:length(sm))-sm(1:length(sm)-1);      %POSITIVE SM INCREMENTS
                 psm(psm<0) = [];
                 dpsoil1=30; % unit mm
                 sm_sum(i,j) = sum(psm).*dpsoil1;% dpsoil1=surface soil moisture depth in NNsm is 3mm.
              end
         end
     end
    prec_sum = nansum(GPM_prec,3); 
     pfrac = sm_sum./prec_sum;
     tauS = -(3/2)./log(pfrac);
    filename = ['/home/xyw/matlabcode/nnsm_layer1/soil_memory/NNsm_tauS_f3_',yrlist{yr},'.mat'];
    save(filename,'tauS','-v7.3');  
    filename = ['/home/xyw/matlabcode/nnsm_layer1/soil_memory/NNsm_pfrac_f3_',yrlist{yr},'.mat'];
    save(filename,'pfrac','-v7.3');  
end
%% imshow
% yrlist=...
%     {'2002','2003','2004','2005','2006',...
%     '2007','2008' ,'2009','2010','2011'...
%      '2012','2013','2014','2015'};
% 
% for yr=1:length(yrlist)
%     load(['/home/xyw/matlabcode/nnsm_layer1/NNsm_tauS_f3_',yrlist{yr},'.mat']);
%     imagesc(tauS),colorbar,caxis([0 0.2]);
%     title(yrlist(yr));
%     pause(1);
% end