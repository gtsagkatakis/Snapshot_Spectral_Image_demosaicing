function [mean_PSNR,mean_SAM]=evaluate_on_SENTINEL(num_band,sz,smp_scenario,num_obs_pxl,GRMR_params)

load('Sentinel_data.mat');

I_HS=I_HS(1:sz(1),1:sz(2),:);
mx=max(max(max(I_HS)));
I_HS=I_HS./mx;
I_HS=I_HS*255;
I_HS=round(I_HS);

[n1,n2,n3]=size(I_HS);

[SMP_seq,FilterPattern_lst]=make_sampling_operators(n1,n2,n3,num_obs_pxl,num_band,smp_scenario,SpectralProfiles);
[I_MOS_seq]=acquire_observations(I_HS,SMP_seq,num_obs_pxl);

SMP_SEQ=SMP_seq;

I_WB_tmp=zeros(sz(1),sz(2),num_band,num_obs_pxl);
I_BTES_tmp=zeros(sz(1),sz(2),num_band,num_obs_pxl);
I_ItSD_tmp=zeros(sz(1),sz(2),num_band,num_obs_pxl);
I_PPID_tmp=zeros(sz(1),sz(2),num_band,num_obs_pxl);


for pp=1:num_obs_pxl
    I_MOS=I_MOS_seq(:,:,pp);
    FilterPattern=cell2mat(FilterPattern_lst(pp));
    
    disp('Running WB');
    I_WB_tmp(:,:,:,pp)=run_WB(I_MOS,FilterPattern,num_band);
    disp('Running BTES');
    I_BTES_tmp(:,:,:,pp)=run_BTES(I_MOS,FilterPattern,num_band,squeeze(I_WB_tmp(:,:,:,pp)));
    disp('Running ItSD');
    I_ItSD_tmp(:,:,:,pp)=run_ItSD(I_MOS,FilterPattern,num_band);
    disp('Running PPID');
    PPI=mean(squeeze(I_WB_tmp(:,:,:,pp)),3);
    I_PPID_tmp(:,:,:,pp)=run_PPID(I_MOS,FilterPattern,num_band,PPI);
end

I_WB=mean(I_WB_tmp,4);
I_BTES=mean(I_BTES_tmp,4);
I_ItSD=mean(I_ItSD_tmp,4);
I_PPID=mean(I_PPID_tmp,4);

% offset=10;
% maxIter=20;
% sgm2=1e2; %weight for distnace estimation
% gamma=0.1; %0.1
% rank_sel=2;

offset=GRMR_params.offset;
maxIter=GRMR_params.maxIter; 
sgm2=GRMR_params.sgm2; 
gamma=GRMR_params.gamma; 
rank_sel=GRMR_params.rank_sel;

disp('Running GRMR');
I_GRMR_rec=run_GRMR_demosaick(I_MOS_seq,SMP_SEQ,num_band,offset,sgm2,maxIter,rank_sel,gamma,I_WB);

for band=1:num_band
    err_GRMR(band)=psnr(squeeze(I_GRMR_rec(:,:,band)),squeeze(I_HS(:,:,band)),256);
    err_PPID(band)=psnr(squeeze(I_PPID(:,:,band)),squeeze(I_HS(:,:,band)),256);
    err_WB(band)=psnr(squeeze(I_WB(:,:,band)),squeeze(I_HS(:,:,band)),256);
    err_ItSD(band)=psnr(squeeze(I_ItSD(:,:,band)),squeeze(I_HS(:,:,band)),256);
    err_BTES(band)=psnr(squeeze(I_BTES(:,:,band)),squeeze(I_HS(:,:,band)),256);
end

for xx=1:size(I_GRMR_rec,1)
    for yy=1:size(I_GRMR_rec,2)
        tmp0=squeeze(I_HS(xx,yy,:))+eps;
        
        tmp1=round(squeeze(I_GRMR_rec(xx,yy,:)))+eps;
        err_SAM_GRMR_temp(xx,yy)=real(hyperSam(tmp1,tmp0));
        
        tmp1=round(squeeze(I_PPID(xx,yy,:)))+eps;
        err_SAM_PPID_temp(xx,yy)=real(hyperSam(tmp1,tmp0));
        
        tmp1=round(squeeze(I_WB(xx,yy,:)))+eps;
        err_SAM_WB_temp(xx,yy)=real(hyperSam(tmp1,tmp0));
        
        tmp1=round(squeeze(I_ItSD(xx,yy,:)))+eps;
        err_SAM_ItSD_temp(xx,yy)=real(hyperSam(tmp1,tmp0));
        
        tmp1=round(squeeze(I_BTES(xx,yy,:)))+eps;
        err_SAM_BTES_temp(xx,yy)=real(hyperSam(tmp1,tmp0));
        
    end
end
err_SAM_GRMR=mean(mean(err_SAM_GRMR_temp));
err_SAM_PPID=mean(mean(err_SAM_PPID_temp));
err_SAM_WB=mean(mean(err_SAM_WB_temp));
err_SAM_ItSD=mean(mean(err_SAM_ItSD_temp));
err_SAM_BTES=mean(mean(err_SAM_BTES_temp));

mean_PSNR=[mean(mean(err_GRMR)),mean(mean(err_PPID)),mean(mean(err_WB)),mean(mean(err_ItSD)),mean(mean(err_BTES))];
mean_SAM=[mean((err_SAM_GRMR)),mean((err_SAM_BTES)),mean((err_SAM_WB)),mean((err_SAM_PPID)),mean((err_SAM_ItSD))];

h1=figure;
h2=figure;
h3=figure;
h4=figure;
h5=figure;

for tt=1:num_band
    figure(h1); imagesc(squeeze(I_GRMR_rec(:,:,tt)),[0,256]); colormap('gray'); title('GRMR');
    figure(h2); imagesc(squeeze(I_BTES(:,:,tt)),[0,256]); colormap('gray');title('BTES');
    figure(h3); imagesc(squeeze(I_WB(:,:,tt)),[0,256]); colormap('gray');title('WB');
    figure(h4); imagesc(squeeze(I_PPID(:,:,tt)),[0,256]); colormap('gray');title('PPID');
    figure(h5); imagesc(squeeze(I_ItSD(:,:,tt)),[0,256]); colormap('gray');title('ItSD');
    
    pause(1);
    
end



