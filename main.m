close all
clearvars
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%exit
%%
num_band=16; %num_band_lst(tt);    %number of bands
sz=[500,500];     %image size
smp_scenario=3; %smp_scenario_lst(qq); %sampling operator, 1=binary, 2=random projections, 3=filter based
num_obs_pxl=2; % num_obs_pxl_lst(zz);  %number of exposures

GRMR_params.offset=5;
GRMR_params.maxIter=5; %was 20
GRMR_params.sgm2=1e1;
GRMR_params.gamma=0.2;
GRMR_params.rank_sel=2;

[mean_PSNR,mean_SAM,exec_time,std_PSNSR,std_SAM]=evaluate_on_CAVE(num_band,sz,smp_scenario,num_obs_pxl,GRMR_params);

fprintf('\n');
fprintf('GMRM: PSNR=%.1f (%.1f), SAM=%.2f (%.2f), Time=%.2f\n',mean_PSNR(1),std_PSNSR(1),mean_SAM(1),std_SAM(1),exec_time(1));
fprintf('BTES: PSNR=%.1f (%.1f), SAM=%.2f (%.2f), Time=%.2f\n',mean_PSNR(2),std_PSNSR(2),mean_SAM(2),std_SAM(2),exec_time(2));
fprintf('WB:   PSNR=%.1f (%.1f), SAM=%.2f (%.2f), Time=%.2f\n',mean_PSNR(3),std_PSNSR(3),mean_SAM(3),std_SAM(3),exec_time(3));
fprintf('PPID: PSNR=%.1f (%.1f), SAM=%.2f (%.2f), Time=%.2f\n',mean_PSNR(4),std_PSNSR(4),mean_SAM(4),std_SAM(4),exec_time(4));
fprintf('ItSD: PSNR=%.1f (%.1f), SAM=%.2f (%.2f), Time=%.2f\n',mean_PSNR(5),std_PSNSR(5),mean_SAM(5),std_SAM(5),exec_time(5));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

num_band=9;    %number of bands
sz=[600,600];     %image size
smp_scenario=3; %sampling operator, 1=binary, 2=random projections, 3=filter based
num_obs_pxl=1;  %number of exposures

GRMR_params.offset=5;
GRMR_params.maxIter=10; %was 20
GRMR_params.sgm2=1e1;
GRMR_params.gamma=0.2;
GRMR_params.rank_sel=2;

[mean_PSNR,mean_SAM]=evaluate_on_SENTINEL(num_band,sz,smp_scenario,num_obs_pxl,GRMR_params);


fprintf('\n');
fprintf('GMRM: PSNR=%.2f, SAM=%.2f\n',mean_PSNR(1),mean_SAM(1));
fprintf('BTES: PSNR=%.2f, SAM=%.2f\n',mean_PSNR(2),mean_SAM(2));
fprintf('WB:   PSNR=%.2f, SAM=%.2f\n',mean_PSNR(3),mean_SAM(3));
fprintf('PPID: PSNR=%.2f, SAM=%.2f\n',mean_PSNR(4),mean_SAM(4));
fprintf('ItSD: PSNR=%.2f, SAM=%.2f\n',mean_PSNR(5),mean_SAM(5));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
disp('IMEC data');
num_band=16;    %number of bands
sz=[1000,1000];     %image size
smp_scenario=3; %sampling operator, 1=binary, 2=random projections, 3=filter based
num_obs_pxl=1;  %number of exposures

GRMR_params.offset=5;
GRMR_params.maxIter=10; %was 20
GRMR_params.sgm2=1e1;
GRMR_params.gamma=0.2;
GRMR_params.rank_sel=2;

RES=evaluate_on_IMEC(num_band,sz,smp_scenario,num_obs_pxl,GRMR_params);
I_HS=RES{1};
I_GRMR_rec=RES{2};
I_BTES=RES{3};
I_WB=RES{4};
I_PPID=RES{5};
I_ItSD=RES{6};
figure;

for tt=1:num_band
    tt
    subplot(2,3,1); imagesc(squeeze(I_HS),[0,256]); colormap('gray'); title('Input Mosaic');
    subplot(2,3,2); imagesc(squeeze(I_GRMR_rec(:,:,tt)),[0,256]); colormap('gray'); title('GRMR');
    subplot(2,3,3); imagesc(squeeze(I_BTES(:,:,tt)),[0,256]); colormap('gray');title('BTES');
    subplot(2,3,4); imagesc(squeeze(I_WB(:,:,tt)),[0,256]); colormap('gray');title('WB');
    subplot(2,3,5); imagesc(squeeze(I_PPID(:,:,tt)),[0,256]); colormap('gray');title('PPID');
    subplot(2,3,6); imagesc(squeeze(I_ItSD(:,:,tt)),[0,256]); colormap('gray');title('ItSD');
    
    pause(1);
    
end

