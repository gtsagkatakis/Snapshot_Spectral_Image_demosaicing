function I_REF=load_hypercube_CAVE(fname,sz,num_band)

for tt=1:num_band
    tmp=double(imread(sprintf('complete_ms_data\\%s_ms\\%s_ms_%02d.png',fname,fname,tt)));
    tmp2=tmp(1:sz(1),1:sz(2));
    I_REF(:,:,tt)=tmp2;
end