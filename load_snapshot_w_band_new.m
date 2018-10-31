function [I1_MOS,I1_SMP_SEQ,band_sorted]=load_snapshot_w_band_new(sz,num_band)

if num_band==25
    load spectral_responses_5x5
%     I1_MOS=double(imread('example1_5x5.pgm'));
    I1_MOS=double(imread('example2_5x5.pgm'));
    CentralWavelengths=round(CentralWavelengths);
elseif num_band==16
    load spectral_responses_4x4
    Wavelength = round(CentralWavelength);
    I1_MOS=double(imread('example1_4x4.pgm'));
%     I1_MOS=double(imread('example2_4x4.pgm'));
    CentralWavelengths=round(CentralWavelength);
end

band_sel=round(CentralWavelengths-400);

S=SpectralProfiles(:,band_sel);

[~,band_sorted]=sort(CentralWavelengths);


for tt=1:num_band
    S(tt,:)=S(tt,:)/sum(S(tt,:));
end


S2=reshape(S,[sqrt(num_band),sqrt(num_band),num_band]);
I1_SMP_SEQ=repmat(S2,[sz(1)/sqrt(num_band),sz(2)/sqrt(num_band),1]);


I1_MOS=I1_MOS(1:sz(1),1:sz(2));
I1_MOS=(I1_MOS./max(max(I1_MOS)))*255;


return
