function I_WB=run_WB(I_MOS,FilterPattern,num_band)

[n1,n2]=size(I_MOS);
I_hat=zeros(n1,n2,num_band);

if num_band==9
   F= (1/9)*[   1,2,3,2,1;...
                2,4,6,4,2;...
                3,6,9,6,3;...
                2,4,6,4,2;...
                1,2,3,2,1];
elseif num_band==16
    F= (1/16)*[1,2,3,4,3,2,1;...
                2,4,6,8,6,4,2;...
                3,6,9,12,9,6,3;...
                4,8,12,16,12,8,4;...
                3,6,9,12,9,6,3;...
                2,4,6,8,6,4,2;...
                1,2,3,4,3,2,1];
elseif  num_band==25
    F= (1/25)*[ 1,  2,  3,  4,  5,  4,  3,  2,  1;...
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


I_3D=zeros(n1,n2,num_band);
for xx=1:n1
    for yy=1:n2
        tmp1=I_MOS(xx,yy);
        tmp2=FilterPattern(xx,yy);
        I_3D(xx,yy,tmp2)=tmp1;
    end
end

I_WB=zeros(n1,n2,num_band);
for tt=1:num_band
    tmp1=squeeze(I_3D(:,:,tt));
    tmp2 = padarray(tmp1,[7,7],'replicate');
    tmp3=conv2(tmp2,F,'same');
    I_WB(:,:,tt)=tmp3(8:end-7,8:end-7);
%     for xx=4:n1-3
%         for yy=4:n2-3
%             
%         end
%     end
end

% for tt=1:num_band
%     tmp=(FilterPattern==tt);
%     I_band=I_MOS.*tmp;
%     tmp2=conv2(I_band,F,'same');
%     I_WB(:,:,tt)=tmp2;
% end


