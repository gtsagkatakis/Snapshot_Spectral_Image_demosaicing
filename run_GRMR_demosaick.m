function I_GMRM_rec=run_GRMR_demosaick(I_MOS_seq,SMP_seq,num_band,offset,sgm2,maxIter,rank_sel,gamma,I_WB)

% step_size=floor(offset);
step_size=floor(offset/2);


[n1,n2,num_obs_pxl]=size(I_MOS_seq);
if num_obs_pxl==1
    [I_MOS_seq,SMP_seq]=replicate_mosaic(I_MOS_seq,SMP_seq,num_obs_pxl,n1,n2,num_band);
    [n1,n2,num_obs_pxl]=size(I_MOS_seq);
end

block_iter=1;
for xx=1:step_size:size(I_MOS_seq,1)-offset+1
    for yy=1:step_size:size(I_MOS_seq,2)-offset+1
        tmp2=I_MOS_seq(xx:xx+offset-1,yy:yy+offset-1,:);
        CurrBlock_lst{block_iter}=reshape(tmp2,[offset*offset,num_obs_pxl]);
        tmp2b=I_WB(xx:xx+offset-1,yy:yy+offset-1,:);
        INT{block_iter}=reshape(tmp2b,[offset*offset,num_band]);
        tmp3=SMP_seq(xx:xx+offset-1,yy:yy+offset-1,:,:);
        tmp4=reshape(tmp3,[offset*offset,num_band,size(tmp3,4)]);
        FilterPattern_lst{block_iter}=tmp4;
        Loc{block_iter}=[xx,yy];
        
        block_iter=block_iter+1;
    end
end

I_rec0=zeros(n1,n2,num_band);

parfor tt=1:numel(FilterPattern_lst)
    L_lst{tt}= make_Laplacian(INT{tt},sgm2);
end

D_rec0=cell(numel(FilterPattern_lst),1);
parfor tt=1:numel(FilterPattern_lst)
    M2=CurrBlock_lst{tt};
    smp2=FilterPattern_lst{tt};
    L=L_lst{tt};
    WB=INT{tt};
    D_rec0{tt}=run_GRMR_final(smp2,M2,L,maxIter,gamma,rank_sel,WB);
end

CNT=zeros(n1,n2,num_band);
for tt=1:numel(FilterPattern_lst)
    tmp=reshape(D_rec0{tt},[offset,offset,num_band]);
    I_rec0(Loc{tt}(1):Loc{tt}(1)+offset-1,Loc{tt}(2):Loc{tt}(2)+offset-1,:)=tmp+I_rec0(Loc{tt}(1):Loc{tt}(1)+offset-1,Loc{tt}(2):Loc{tt}(2)+offset-1,:);
    CNT(Loc{tt}(1):Loc{tt}(1)+offset-1,Loc{tt}(2):Loc{tt}(2)+offset-1,:)=CNT(Loc{tt}(1):Loc{tt}(1)+offset-1,Loc{tt}(2):Loc{tt}(2)+offset-1,:)+ones(offset,offset,num_band);
end
I_GMRM_rec=I_rec0./CNT;

I_GMRM_rec(1:step_size,:,:)=I_WB(1:step_size,:,:);
I_GMRM_rec(:,1:step_size,:)=I_WB(:,1:step_size,:);
I_GMRM_rec(end-step_size-1:end,:,:)=I_WB(end-step_size-1:end,:,:);
I_GMRM_rec(:,end-step_size-1:end,:)=I_WB(:,end-step_size-1:end,:);
