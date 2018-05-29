close all
clearvars
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
addpath(genpath('\\complete_ms_data'));

num_band=16;    %number of bands
sz=[40,40];     %image size
smp_scenario=3; %sampling operator, 1=binary, 2=random projections, 3=filter based
num_obs_pxl=2;  %number of exposures

GRMR_params.offset=10;
GRMR_params.maxIter=20; %was 20
GRMR_params.sgm2=1e2; 
GRMR_params.gamma=0.1; 
GRMR_params.rank_sel=2;

[mean_PSNR,mean_SAM]=evaluate_on_CAVE(num_band,sz,smp_scenario,num_obs_pxl,GRMR_params);

fprintf('\n');
fprintf('GMRM: PSNR=%.2f, SAM=%.2f\n',mean_PSNR(1),mean_SAM(1));
fprintf('BTES: PSNR=%.2f, SAM=%.2f\n',mean_PSNR(2),mean_SAM(2));
fprintf('WB:   PSNR=%.2f, SAM=%.2f\n',mean_PSNR(3),mean_SAM(3));
fprintf('PPID: PSNR=%.2f, SAM=%.2f\n',mean_PSNR(4),mean_SAM(4));
fprintf('ItSD: PSNR=%.2f, SAM=%.2f\n',mean_PSNR(5),mean_SAM(5));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

num_band=9;    %number of bands
sz=[90,90];     %image size
smp_scenario=3; %sampling operator, 1=binary, 2=random projections, 3=filter based
num_obs_pxl=2;  %number of exposures

GRMR_params.offset=10;
GRMR_params.maxIter=20; %was 20
GRMR_params.sgm2=1e2; 
GRMR_params.gamma=0.1; 
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
num_band=25;    %number of bands
sz=[200,200];     %image size
smp_scenario=3; %sampling operator, 1=binary, 2=random projections, 3=filter based
num_obs_pxl=1;  %number of exposures

GRMR_params.offset=10;
GRMR_params.maxIter=20; %was 20
GRMR_params.sgm2=1e2; 
GRMR_params.gamma=0.1; 
GRMR_params.rank_sel=2;

evaluate_on_IMEC(num_band,sz,smp_scenario,num_obs_pxl,GRMR_params);
