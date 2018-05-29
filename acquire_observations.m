function [I_smp]=acquire_observations(I_HS,SMP_seq,num_obs_pxl)

[n1,n2,~]=size(I_HS);

for xx=1:n1
    for yy=1:n2
        tmp1=squeeze(I_HS(xx,yy,:));
        for tt=1:num_obs_pxl
           tmp2= squeeze(SMP_seq(xx,yy,:,tt));
           I_smp(xx,yy,tt)=tmp1'*tmp2; 
        end
    end
end

       