%% Script to estimate tauL and tauS from MTDCA SSM data (2015-2020)
 % Author: Qinghe
 
 %tauL
 load('/mnt/WCL-S920B_home/qinghe/reanaly_compar/latlonP.mat');
 
%  yrlist={'2002','2003','2004','2005','2006'};
%  yrlist={'2007','2008','2009','2010','2011'};
%  yrlist={'2012','2013','2014','2015','2016'};
yrlist={'2022','2023'};

 for yr = 1:length(yrlist)
     
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

     doy = 1:size(SSM_Data_us,3);
     
     num_drydown_all =nan(size(lat,1),size(lat,2));
     sm_prior_all    =nan(size(lat,1),size(lat,2));
     lamuda_mean_all =nan(size(lat,1),size(lat,2));
     
     psm    = nan(1000000,5);
     lamuda = psm;
     theta0 = psm;  % KAM edit: include theta0 (different to psm)
     thetaf = psm;  % KAM edit: include thetaf (lower bound parameter)
     
     R2       = nan(1000000,3);
     dd_dates = nan(1000000,4);
     
     ii = 1;
     for i = 1:size(lat,1)
%          disp(['MTDCA, ',yrlist{yr},': ',num2str(i)])
        disp(['NNsm_year, ',yrlist{yr},': ',num2str(i)])
         for j = 1:size(lat,2)
             sm = reshape(SSM_Data_us(i,j,:),1,size(SSM_Data_us,3));
             sm(sm<0.04) = nan;
             sm(sm>0.7) = nan;
             t_threshold = (max(sm)-min(sm))*0.1;
             if ~isempty(find(sm>0, 1))
                 % info: start doy (col 1), end doy (col 2), num of eff
                 % observations (col 3), dry down days (col 4)
                 [info,vol]=detect_drydownn(sm, doy, t_threshold);
                 
                 [n,m]=size(info);
                 if n>0
                    slp_dd=zeros(n,1);
                    slp_dd(:)=nan;
                    prior_mm2=slp_dd;
                    prior_mm1=slp_dd;
                    lengthofdd=prior_mm2;
                    
                    for id=1:n
                        d_st=find(doy==info(id,1));%dry down start date
                        d_ed=find(doy==info(id,2));%dry down end date
                        x=doy(d_st:d_ed);%date
                        y=sm(d_st:d_ed);%sm
                        A=find(y>0);%remove nan data
                        x=x(A);
                        x=x-x(1)+1;
                        x= reshape(x,length(x),1);
                        y=double(y(A)');
                        nn=sum(~isnan(x));
                        y0=sm(d_st);
                        
                        if nn>3%select eff obs greater than 3
                            if nn<11
                                % KAM edit: add an additional parameter to
                                % the model to be fit
                                p=fittype('y0*exp(-x*c)+yf','independent','x');
                                %p=fittype('y0*exp(-x*c)','independent','x');
                                % KAM edit: also output goodness-of-fit
                                % statistics; set bounds on fitting yf
                                %[f,gof]=fit(x,y,p,'StartPoint', [0.0001, sm(d_st)]);
                                miny = nanmin(y);
                                if  miny-nanmin(sm) <1e-05
                                     miny = nanmin(sm)+1e-04;
                                end
                                [f,gof]=fit(x,y,p,'StartPoint', [0.0001, sm(d_st), nanmin(y)],...
                                    'Lower',[0, 0, nanmin(sm)],'Upper',[1e3, 1,  miny]);
                                % KAM: Coefficient order is c, y0, yf
                                % KAM: consistent with Shellito et al.
                                % (2016), we constrain yf below the lowest
                                % SM value observed in the individual
                                % drydown, but above the minimum SM value
                                % observed over the entire time period.
                                ci = confint(f); %QING:confidence intervals
                                %f=fit(x,y,p,'StartPoint', [0.0001, sm(d_st)]);
                                lengthofdd(id)=nn;
                                slp_dd(id)=f.c;
                                prior_mm1(id)=sm(d_st);
                                lamuda(ii,3)=f.c;
                                lamuda(ii,4)=ci(1,1);
                                lamuda(ii,5)=ci(2,1);
                                psm(ii,3)=sm(d_st);
                                theta0(ii,3) = f.y0;
                                theta0(ii,4) = ci(1,2);
                                theta0(ii,5) = ci(2,2);
                                thetaf(ii,3) = f.yf;
                                thetaf(ii,4) = ci(1,3);
                                thetaf(ii,5) = ci(2,3);
                                lamuda(ii,1)=i;
                                lamuda(ii,2)=j;
                                psm(ii,1)=i;
                                psm(ii,2)=j;
                                theta0(ii,1)=i;
                                theta0(ii,2)=j;
                                thetaf(ii,1)=i;
                                thetaf(ii,2)=j;
                                R2(ii,1) = i;
                                R2(ii,2) = j;
                                R2(ii,3) = gof.rsquare;
                                dd_dates(ii,1) = i;
                                dd_dates(ii,2) = j;
                                dd_dates(ii,3) = d_st;
                                dd_dates(ii,4) = d_ed;
                                clear p f
                                ii=ii+1;
                            end
                        end
                    end
                    
                    num_drydown_all(i,j)=sum(~isnan(slp_dd));
                    sm_prior_all(i,j)=nanmean(prior_mm1);
                    lamuda_mean_all(i,j)=nanmean(slp_dd);
                    
                 end
             end
         end
     end
     
     filename = ['/home/xyw/matlabcode/nnsm_layer1/NNsm_tauL_f3_',yrlist{yr},'.mat'];
     save(filename,'sm_prior_all', 'lamuda_mean_all', ...
                   'num_drydown_all','lamuda',...
                   'psm','R2','theta0',...
                   'thetaf','dd_dates','-v7.3');           
 end
      

     