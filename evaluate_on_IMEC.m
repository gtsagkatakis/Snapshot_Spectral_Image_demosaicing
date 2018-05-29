function evaluate_on_IMEC(num_band,sz,smp_scenario,num_obs_pxl,GRMR_params)


[I_HS,I1_SMP_SEQ]=load_snapshot_w_band_new(sz,num_band);


[I_MOS_seq]=simulate_video(I_HS,I1_SMP_SEQ,num_obs_pxl,num_band);


if num_band==25
    load('spectral_responses_5x5.mat');
elseif num_band==16
    load('spectral_responses_4x4.mat');
    CentralWavelengths=CentralWavelength;
else
    disp('Error');
end

temp2=sort( round(CentralWavelengths))-400;
SpectralProfiles=SpectralProfiles(:,temp2);
SpectralProfiles=rot90(SpectralProfiles);


mx=max(max(max(I_HS)));
I_HS=I_HS./mx;
I_HS=I_HS*255;
I_HS=round(I_HS);

[n1,n2,n3]=size(I_HS);

[SMP_seq,FilterPattern_lst]=make_sampling_operators(n1,n2,n3,num_obs_pxl,num_band,smp_scenario,SpectralProfiles);

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
    disp('Running PPID');
    PPI=mean(squeeze(I_WB_tmp(:,:,:,pp)),3);
    I_PPID_tmp(:,:,:,pp)=run_PPID(I_MOS,FilterPattern,num_band,PPI);
    disp('Running BTES');
    I_BTES_tmp(:,:,:,pp)=run_BTES(I_MOS,FilterPattern,num_band,squeeze(I_WB_tmp(:,:,:,pp)));
    disp('Running ItSD');
    I_ItSD_tmp(:,:,:,pp)=run_ItSD(I_MOS,FilterPattern,num_band);
    
end

I_WB=mean(I_WB_tmp,4);
I_BTES=mean(I_BTES_tmp,4);
I_ItSD=mean(I_ItSD_tmp,4);
I_PPID=mean(I_PPID_tmp,4);

% offset=10;
% maxIter=20; %was 20
% sgm2=1e2; 
% gamma=0.1; 
% rank_sel=2;

offset=GRMR_params.offset;
maxIter=GRMR_params.maxIter; 
sgm2=GRMR_params.sgm2; 
gamma=GRMR_params.gamma; 
rank_sel=GRMR_params.rank_sel;


disp('Running GRMR');
I_GRMR_rec=run_GRMR_demosaick(I_MOS_seq,SMP_SEQ,num_band,offset,sgm2,maxIter,rank_sel,gamma,I_WB);



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



