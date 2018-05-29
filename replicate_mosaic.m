function [I_MOS_seq2,SMP_seq2]=replicate_mosaic(I_MOS_seq,SMP_seq,num_obs_pxl,n1,n2,num_band)

if num_obs_pxl==1
    
    I_MOS_seq2=repmat(I_MOS_seq,[1,1,9]);
    SMP_seq2=repmat(SMP_seq,[1,1,1,9]);
    idx_x=2;
    
    for xx=2:1:n1-1
        idx_y=2;
        for yy=2:1:n2-1
            tmp=I_MOS_seq(xx-1:xx+1,yy-1:yy+1);
            I_MOS_seq2(idx_x,idx_y,1:9)=tmp(:);
            tmp2=SMP_seq(xx-1:xx+1,yy-1:yy+1,:);
            SMP_seq2(idx_x,idx_y,:,:)=reshape(tmp2,[9,num_band])';
            idx_y=idx_y+1;
        end
        idx_x=idx_x+1;
        
    end
end
