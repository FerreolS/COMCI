Imphase = (imread('Data/Simus/GrndTruth/OrionNebula_Amp_Ph.tif','Index',2));
ImMod = (imread('Data/Simus/GrndTruth/OrionNebula_Amp_Ph.tif','Index',1));

%w1 = window(@tukeywin,18000,0.002);
w1 =  1;

truth = (1.0-1./255.*double((double(ImMod).*w1) .* w1')) .* exp(1i.* pi.* double(Imphase)./180);
clear ImMod Imphase
padsz=1080;
truth = padarray(truth,[padsz, padsz],'symmetric');
szt = size(truth);