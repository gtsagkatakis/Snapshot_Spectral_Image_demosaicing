function I_REC=ItSD(I_MOS,FilterPattern,num_band)

I_int=run_WB(I_MOS,FilterPattern,num_band);

H=(1/16)*[ 1, 2, 3, 4, 3, 2, 1;...
    2, 4, 6, 8, 6, 4, 2;...
    3, 6, 9, 12,9, 6, 3;...
    4, 8, 12,16,12,8, 4;...
    3, 6, 9, 12,9, 6, 3;...
    2, 4, 6, 8, 6, 4, 2;...
    1, 2, 3, 4, 3, 2, 1];

Imean=I_int;
for iter=1:5
    Imean=mean(Imean,3);
    
    for tt=1:num_band
        smp_tmp=double(FilterPattern==tt);
        I_k=I_MOS.*smp_tmp;
        I_d=squeeze(I_int(:,:,tt)).*smp_tmp;
        D_k=I_k-I_d;
        D_k_hat=conv2(D_k,H,'same');
        I_REC(:,:,tt)=Imean+D_k_hat;
    end
    
end