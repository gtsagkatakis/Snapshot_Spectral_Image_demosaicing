function [SMP_seq,FilterPattern_lst]=make_sampling_operators(n1,n2,n3,num_obs_pxl,num_band,smp_scenario,SpectralProfiles)

for tt=1:num_band
    
    tmp=(rand(1,num_band));
    tmp=tmp/sum(tmp);
    rand_pattern(tt,:)=tmp;
end


for ww=1:num_obs_pxl
    ptn1=randperm(num_band);
    ptn2=reshape(ptn1,[sqrt(num_band),sqrt(num_band)]);
    r=n1/sqrt(num_band);
    FilterPattern=repmat(ptn2,[r,r]);
    FilterPattern_lst{ww}=FilterPattern;
end



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

