function [I_final]=run_PPID(I_MOS,FilterPattern,num_band,PPI)

[n1,n2]=size(I_MOS);
% M=(1/64).*[ 1, 2, 2, 2, 1;...
%     2, 4, 4, 4, 2;...
%     2, 4, 4, 4, 2;...
%     2, 4, 4, 4, 2;...
%     1, 2, 2, 2, 1];
% 
% I_hat=conv2(I_MOS,M,'same');
% 
 
I_hat=PPI; 
 


I_D=zeros(n1,n2,num_band);
for tt=1:num_band
    tmp=(FilterPattern==tt);
    I_D(:,:,tt)=(I_MOS-PPI).*tmp;
    
end
 
for xx=1:n1
    for yy=1:n2
        lst=[-4,-4;-4,0;-4,4;0,4;4,4;4,0; 0,-4;4,-4];
        Neigh_coords=repmat([xx,yy],[8,1])+lst;
        Neigh_coords(:,2)=max(1,Neigh_coords(:,2));
        Neigh_coords(:,2)=min(n2,Neigh_coords(:,2));
        Neigh_coords(:,1)=max(1,Neigh_coords(:,1));
        Neigh_coords(:,1)=min(n1,Neigh_coords(:,1));

        for qq=1:8
            Neigh_val=I_hat(Neigh_coords(qq,1),Neigh_coords(qq,2));
            gm(qq)=1/(1+abs(I_hat(xx,yy)-Neigh_val));

        end
        gm_lst{xx,yy}=gm;
    end

end

if num_band==9
   Q= (1/9)*[   1,2,3,2,1;...
                2,4,6,4,2;...
                3,6,9,6,3;...
                2,4,6,4,2;...
                1,2,3,2,1];
elseif num_band==16
    Q= (1/16)*[1,2,3,4,3,2,1;...
                2,4,6,8,6,4,2;...
                3,6,9,12,9,6,3;...
                4,8,12,16,12,8,4;...
                3,6,9,12,9,6,3;...
                2,4,6,8,6,4,2;...
                1,2,3,4,3,2,1];
elseif  num_band==25
    Q= (1/25)*[ 1,  2,  3,  4,  5,  4,  3,  2,  1;...
                2,  4,  6,  8, 10,  8,  6,  4,  2;...
                3,  6,  9, 12, 15, 12,  9,  6,  3;...
                4,  8, 12, 16, 20, 16, 12,  8,  4;...
                5,  10,15, 20, 25, 20, 15,  10, 5;...
                4,  8, 12, 16, 20, 16, 12,  8,  4;...
                3,  6,  9, 12, 15, 12,  9,  6,  3;...
                2,  4,  6,  8, 10,  8,  6,  4,  2;...
                1,  2,  3,  4,  5,  4,  3,  2,  1];
    
else
   return 
end



for xx=1:n1
    for yy=1:n2
        gm=gm_lst{xx,yy};
        Gp =    [gm(1)*eye(3) repmat(gm(2),[3,1]) gm(3)*eye(3);...
                 repmat(gm(8),[1,3]) 1  repmat(gm(4),[1,3]);...
                 gm(7)*eye(3) repmat(gm(6),[3,1]) gm(5)*eye(3)];
        for qq=1:sqrt(num_band)
            for zz=1:sqrt(num_band)
                H_loc(qq,zz)=Q(qq,zz)*Gp(qq,zz);
            end
        end
        H{xx,yy}=H_loc;
    end
end

for tt=1:num_band
    I_band=squeeze(I_D(:,:,tt));
    tmp2=conv2(I_band,Q,'same');
    I_D2(:,:,tt)=tmp2;
    I_final(:,:,tt)=PPI+I_D2(:,:,tt);
end

% figure;
% for tt=1:16
    
    %     I_final(:,:,tt)=I_final(:,:,tt)*mx;
%     imagesc(squeeze(I_final(:,:,tt)));
%     img1=squeeze(I_final(:,:,tt));
%     img2=squeeze(I_HS(:,:,tt));
%     error_band(tt)=psnr(img1,img2,255);
%     title(sprintf('Err %f',error_band(tt)));
%     pause(0.1);
% end
% mean(error_band)



