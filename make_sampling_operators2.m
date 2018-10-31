function [SMP_seq,FilterPattern_lst]=make_sampling_operators2(n1,n2,n3,num_obs_pxl,num_band,smp_scenario,SpectralProfiles)

for tt=1:num_band
    tmp=(rand(1,num_band));
    tmp=tmp/sum(tmp);
    rand_pattern(tt,:)=tmp;
end


if num_band==9
    ptn1=[4     3     6     9     1     8     2     5     7];
elseif num_band==16
    ptn1=[12    16     1    13    15     9     4     6    11    14     5     7     2     8     3    10];
elseif num_band==25
     ptn1=[9 4 25 2 13 14 12 3 5 18 1 22 11 20 23 7 8 10 6 19 24 21 15 17 16];
else
    ptn1=randperm(num_band) ;
end
ptn2=reshape(ptn1,[sqrt(num_band),sqrt(num_band)]);

for ww=1:num_obs_pxl
    if mod(ww,sqrt(n3))==1
        ptn2=circshift(ptn2,1,1);  
    else
        ptn2=circshift(ptn2,1,2);
    end
    r=n1/sqrt(num_band);
    FilterPattern=repmat(ptn2,[r,r]);
    FilterPattern_lst{ww}=FilterPattern;
end

% for tt=1:25
%     tmp(tt)=FilterPattern_lst{tt}(1,1);
% end


SMP_seq=zeros(n1,n2,num_band,num_obs_pxl);

for xx=1:n1
    for yy=1:n2
        for qq=1:num_obs_pxl
            FilterPattern=FilterPattern_lst{qq};
            selected_band=FilterPattern(xx,yy);
            
            tmp=zeros(1,n3);
            switch smp_scenario
                case 1
                    tmp(selected_band)=1;
%                     tmp=conv(tmp,gausswin(num_band,20),'same');
%                     tmp=tmp./sum(tmp);
                case 2
                    tmp=rand_pattern(selected_band,:);
                case 3
                    selected_profile=SpectralProfiles(selected_band,:);
                    selected_profile=selected_profile./sum(selected_profile);
                    tmp=selected_profile;
            end
            
            SMP_seq(xx,yy,:,qq)=tmp;
        end
    end
end

