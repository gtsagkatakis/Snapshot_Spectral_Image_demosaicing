function [mean_PSNR,mean_SAM, exec_time,std_PSNSR,std_SAM]=evaluate_on_CAVE(num_band,sz,smp_scenario,num_obs_pxl,GRMR_params)

fnames=dir('complete_ms_data');

for zz=1:numel(fnames)-2
    fname=fnames(zz+2).name;
    
    fprintf('%s, %d out of %d\n',fname,zz,numel(fnames)-2);
    
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
    
    
    [SMP_seq,FilterPattern_lst]=make_sampling_operators2(n1,n2,n3,num_obs_pxl,num_band,smp_scenario,SpectralProfiles);
    [I_MOS_seq]=acquire_observations(I_HS,SMP_seq,num_obs_pxl);
    
    SMP_SEQ=SMP_seq;
    
    I_WB_tmp=zeros(sz(1),sz(2),num_band,num_obs_pxl);
    I_BTES_tmp=zeros(sz(1),sz(2),num_band,num_obs_pxl);
    I_ItSD_tmp=zeros(sz(1),sz(2),num_band,num_obs_pxl);
    I_PPID_tmp=zeros(sz(1),sz(2),num_band,num_obs_pxl);
    
    disp('Running WB');
    tic;
    for pp=1:num_obs_pxl
        I_MOS=I_MOS_seq(:,:,pp);
        FilterPattern=cell2mat(FilterPattern_lst(pp));
        I_WB_tmp(:,:,:,pp)=run_WB(I_MOS,FilterPattern,num_band);
    end
    I_WB=mean(I_WB_tmp,4);
    WB_toc=toc;
    
    offset=GRMR_params.offset;
    maxIter=GRMR_params.maxIter;
    sgm2=GRMR_params.sgm2;
    gamma=GRMR_params.gamma;
    rank_sel=GRMR_params.rank_sel;
    disp('Running GRMR');
    tic;
    I_GRMR_rec=run_GRMR_demosaick(I_MOS_seq,SMP_SEQ,num_band,offset,sgm2,maxIter,rank_sel,gamma,I_WB);
    GRMR_toc=toc;
    
    
    disp('Running PPID');
    tic
    for pp=1:num_obs_pxl
        I_MOS=I_MOS_seq(:,:,pp);
        FilterPattern=cell2mat(FilterPattern_lst(pp));
        PPI=mean(squeeze(I_WB_tmp(:,:,:,pp)),3);
        I_PPID_tmp(:,:,:,pp)=run_PPID(I_MOS,FilterPattern,num_band,PPI);
    end
    I_PPID=mean(I_PPID_tmp,4);
    PPID_toc=toc;
    
    disp('Running BTES');
    tic;
    for pp=1:num_obs_pxl
        I_MOS=I_MOS_seq(:,:,pp);
        FilterPattern=cell2mat(FilterPattern_lst(pp));
        I_BTES_tmp(:,:,:,pp)=run_BTES(I_MOS,FilterPattern,num_band,squeeze(I_WB_tmp(:,:,:,pp)));
    end
    I_BTES=mean(I_BTES_tmp,4);
    BTES_toc=toc;
    I_BTES=I_WB;
    BTES_toc=WB_toc;
    
    disp('Running ItSD');
    tic;
    for pp=1:num_obs_pxl
        I_MOS=I_MOS_seq(:,:,pp);
        FilterPattern=cell2mat(FilterPattern_lst(pp));
        I_ItSD_tmp(:,:,:,pp)=run_ItSD(I_MOS,FilterPattern,num_band);
    end
    I_ItSD=mean(I_ItSD_tmp,4);
    ItSD_toc=toc;
    %
    
    
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
    
    
    h1=figure;
    h2=figure;
    h3=figure;
    h4=figure;
    h5=figure;
    h6=figure;
    for tt=1:9
        figure(h1);  imagesc(squeeze(I_HS(:,:,tt)),[0,255]); title('TRUE'); %colormap('gray');
        figure(h2);  imagesc(squeeze(I_GRMR_rec(:,:,tt)),[0,255]); title('GRMR'); %colormap('gray');
        figure(h3);  imagesc(squeeze(I_BTES(:,:,tt)),[0,255]); title('BTES'); %colormap('gray');
        figure(h4);  imagesc(squeeze(I_WB(:,:,tt)),[0,255]); title('WB'); %colormap('gray');
        figure(h5);  imagesc(squeeze(I_PPID(:,:,tt)),[0,255]); title('PPID'); %colormap('gray');
        figure(h6);  imagesc(squeeze(I_ItSD(:,:,tt)),[0,255]); title('ItSD'); %colormap('gray');
        
        pause(1);
    end
    
end

exec_time=[GRMR_toc,BTES_toc,WB_toc,PPID_toc,ItSD_toc];
mean_PSNR=[mean(mean(err_GRMR)),mean(mean(err_BTES)),mean(mean(err_WB)),mean(mean(err_PPID)),mean(mean(err_ItSD))];
mean_SAM=[mean((err_SAM_GRMR)),mean((err_SAM_BTES)),mean((err_SAM_WB)),mean((err_SAM_PPID)),mean((err_SAM_ItSD))];
std_PSNSR=[std(mean(err_GRMR)),std(mean(err_BTES)),std(mean(err_WB)),std(mean(err_PPID)),std(mean(err_ItSD))];
std_SAM=[std((err_SAM_GRMR)),std((err_SAM_BTES)),std((err_SAM_WB)),std((err_SAM_PPID)),std((err_SAM_ItSD))];
