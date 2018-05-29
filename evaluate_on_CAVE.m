function [mean_PSNR,mean_SAM]=evaluate_on_CAVE(num_band,sz,smp_scenario,num_obs_pxl,GRMR_params)

for zz=1:4
    
    if zz==1
        fname='flowers';
    elseif zz==2
        fname='beads';
    elseif zz==3
        fname='chart_and_stuffed_toy';
    elseif zz==4
        fname='stuffed_toys';
    else
        disp('Error in file');
        return
    end
    
    I_HS=load_hypercube_CAVE(fname,sz,num_band);
    
    mx=max(max(max(I_HS)));
    I_HS=I_HS./mx;
    I_HS=I_HS*255;
    I_HS=round(I_HS);
    
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
    
    
    offset=GRMR_params.offset;
    maxIter=GRMR_params.maxIter;
    sgm2=GRMR_params.sgm2;
    gamma=GRMR_params.gamma;
    rank_sel=GRMR_params.rank_sel;
    disp('Running GRMR');
    I_GRMR_rec=run_GRMR_demosaick(I_MOS_seq,SMP_SEQ,num_band,offset,sgm2,maxIter,rank_sel,gamma,I_WB);
    
    for band=1:num_band
        err_GRMR(zz,band)=psnr(squeeze(I_GRMR_rec(:,:,band)),squeeze(I_HS(:,:,band)),256);
        err_PPID(zz,band)=psnr(squeeze(I_PPID(:,:,band)),squeeze(I_HS(:,:,band)),256);
        err_WB(zz,band)=psnr(squeeze(I_WB(:,:,band)),squeeze(I_HS(:,:,band)),256);
        err_ItSD(zz,band)=psnr(squeeze(I_ItSD(:,:,band)),squeeze(I_HS(:,:,band)),256);
        err_BTES(zz,band)=psnr(squeeze(I_BTES(:,:,band)),squeeze(I_HS(:,:,band)),256);
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
    err_SAM_GRMR(zz)=mean(mean(err_SAM_GRMR_temp));
    err_SAM_PPID(zz)=mean(mean(err_SAM_PPID_temp));
    err_SAM_WB(zz)=mean(mean(err_SAM_WB_temp));
    err_SAM_ItSD(zz)=mean(mean(err_SAM_ItSD_temp));
    err_SAM_BTES(zz)=mean(mean(err_SAM_BTES_temp));
end

mean_PSNR=[mean(mean(err_GRMR)),mean(mean(err_BTES)),mean(mean(err_WB)),mean(mean(err_PPID)),mean(mean(err_ItSD))];
mean_SAM=[mean((err_SAM_GRMR)),mean((err_SAM_BTES)),mean((err_SAM_WB)),mean((err_SAM_PPID)),mean((err_SAM_ItSD))];
