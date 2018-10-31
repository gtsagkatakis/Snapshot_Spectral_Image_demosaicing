function [I_final]=run_BTES(I_MOS,FilterPattern,num_band,I_HS)

[n1,n2]=size(I_MOS);
I_3D=I_HS;
for xx=1:n1
    for yy=1:n2
        tmp1=I_MOS(xx,yy);
        tmp2=FilterPattern(xx,yy);
        I_3D(xx,yy,tmp2)=tmp1;
    end
end

I_final=I_3D;

for qq=1:num_band
    for tt=1:4
        switch tt
            case 1
                S{qq,tt}=(FilterPattern==qq);
                Y = circshift(S{qq,tt},2,1);
                Y = circshift(Y,2,2);
                Snew{qq,tt}=Y;
            case 2
                S{qq,tt}=S{qq,tt-1}+Snew{qq,tt-1};
                Y = circshift(S{qq,tt},2,1);
                Snew{qq,tt} = Y;
            case 3
                S{qq,tt}=S{qq,tt-1}+Snew{qq,tt-1};
                Y = circshift(S{qq,tt},1,1);
                Y = circshift(Y,1,2);
                Snew{qq,tt} = Y;
            case 4
                S{qq,tt}=S{qq,tt-1}+Snew{qq,tt-1};
                Y = circshift(S{qq,tt},1,1);
                Snew{qq,tt} = Y;
        end
    end
end


for tt=1:num_band
    I_band=squeeze(I_3D(:,:,tt));
    for iter=1:4
        S1=Snew{tt,iter};
        for xx=5:n1-4
            for yy=5:n2-4
                if S1(xx,yy)>0
                    x_lst=xx-2:xx+2;
                    y_lst=yy-2:yy+2;
                    tmp1457=I_band(x_lst,y_lst);
                    a_lst=zeros(5,5);
                    val_lst=zeros(5,5);
                    idx=1;
                    for ff=1:numel(x_lst)
                        for dd=1:numel(y_lst)
                            if I_band(x_lst(ff),y_lst(dd))~=0 && isempty(find(isnan(I_band(x_lst(ff),y_lst(dd)))))% && x_lst(ff)~=xx && y_lst(dd)~=yy
                                tmp0=I_band(x_lst(ff),y_lst(dd));
                                tmp1=tmp0;
                                tmp2=tmp0;
                                tmp3=tmp0;
                                tmp4=tmp0;
                                tmp2454=1/(1+abs(tmp1-tmp0)+abs(tmp2-tmp0)+abs(tmp3-tmp0)+abs(tmp4-tmp0));
                                if ~isempty(find(isnan(tmp2454)))
                                    disp('a')
                                end
                                
                                if ~isempty(find(isnan(tmp0)))
                                    disp('b')
                                end
                                
                                a_lst(ff,dd)=tmp2454;
                                val_lst(ff,dd)=tmp0;
                                idx=idx+1;
                                if ((isnan( a_lst(ff,dd))))
                                    disp('c')
                                end
                                
                            end
                        end
                    end
                    tmp35456=sum(a_lst(:));
                    a_lst2=a_lst./(tmp35456+eps);
                    I_band(xx,yy)=round(a_lst2(:)'*val_lst(:));
                    if ((isnan(I_band(xx,yy))))
                        disp('c')
                    end
                    
                    
                end
            end
        end
    end
    I_final(:,:,tt)=I_band;
end



