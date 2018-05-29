function [I_MOS_seq3, SMP_seq_Full3]=simulate_video(I_MOS_1,SMP_seq_Full,num_frames,num_band)

[n1,n2]=size(I_MOS_1);
step_size=sqrt(num_frames);

I_MOS_seq3=zeros(n1/step_size,n2/step_size,num_frames);
SMP_seq_Full3=zeros(n1/step_size,n2/step_size,num_band,num_frames);

idx_xx=1;
idx_yy=1;
for xx=1:step_size:n1-step_size+1
    for yy=1:step_size:n2-step_size+1
        if num_frames==1
            tmp1=I_MOS_1(xx,yy);
            I_MOS_seq3(idx_xx,idx_yy,:)=tmp1(:);
            tmp2=SMP_seq_Full(xx:xx,yy:yy,:);
            SMP_seq_Full3(idx_xx,idx_yy,:,:)=reshape(tmp2,1,num_band)';
        end
        if num_frames==4
            tmp1=I_MOS_1(xx:xx+1,yy:yy+1);
            I_MOS_seq3(idx_xx,idx_yy,:)=tmp1(:);
            tmp2=SMP_seq_Full(xx:xx+1,yy:yy+1,:);
            SMP_seq_Full3(idx_xx,idx_yy,:,:)=reshape(tmp2,4,num_band)';
        end
        
        if num_frames==9
            tmp1=I_MOS_1(xx:xx+2,yy:yy+2);
            I_MOS_seq3(idx_xx,idx_yy,:)=tmp1(:);
            tmp2=SMP_seq_Full(xx:xx+2,yy:yy+2,:);
            SMP_seq_Full3(idx_xx,idx_yy,:,:)=reshape(tmp2,9,num_band)';
        end
        
         if num_frames==16
            tmp1=I_MOS_1(xx:xx+3,yy:yy+3);
            I_MOS_seq3(idx_xx,idx_yy,:)=tmp1(:);
            tmp2=SMP_seq_Full(xx:xx+3,yy:yy+3,:);
            SMP_seq_Full3(idx_xx,idx_yy,:,:)=reshape(tmp2,16,num_band)';
         end
        
        if num_frames==25
            tmp1=I_MOS_1(xx:xx+4,yy:yy+4);
            I_MOS_seq3(idx_xx,idx_yy,1:25)=tmp1(:);
            tmp2=SMP_seq_Full(xx:xx+4,yy:yy+4,:);
            SMP_seq_Full3(idx_xx,idx_yy,:,:)=reshape(tmp2,25,num_band)';
            
        end
        idx_yy=idx_yy+1;
    end
    idx_yy=1;
    idx_xx=idx_xx+1;
end