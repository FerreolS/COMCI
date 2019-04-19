szt = size(truth);

nAngle = 1;
Iangle = zeros(2,nAngle);
lambda= 530e-9; % wavelenght
n0 = 1.0;
z = 250*1e-6;
k = n0*2.*pi/lambda;
%dxy = 5.2e-6; % Camera  pixel size
dxy= 1.04e-6 / 10;

w = padarray( window(@tukeywin,szt(1)/2,0.), szt(1)/4,0.);
truth = (truth.*w).*w';

Frwrd = LinOpPropagator(szt,lambda, n0, z,dxy,  Iangle, 1,1,'BLAS2');
truePropagatedField = Frwrd* (truth);

rngt = ((szt(1) - 2520)/2+1): (szt(1) + 2520)/2;
trueIntensity = abs(truePropagatedField(rngt,rngt)).^2;
%%
clear truePropagatedField;
sr=[1:7 9 10 12 14 15 18 20];  %downsampling factor
SNR = [0 0.001 0.002 0.005 0.0075 0.01 0.015 0.02 0.05 0.075 0.1 0.2 0.5 1] ;
wdSz = [0 0.01 0.02 0.05 0.1 0.15 0.2 0.3 0.4 0.5 0.75 1. 1.25 1.5  1.75 2 3 4 5 6 7]*1000;

PSNR = zeros(numel(wdSz),numel(sr),numel(SNR));
PSNR2 = zeros(numel(wdSz),numel(sr),numel(SNR));
PSNR3 = zeros(numel(wdSz),numel(sr));


%%
for n=1:numel(wdSz)
    %%  w1 = padarray( window(@tukeywin,1000*n,1000*n/400), (20000 - 1000*n)/2.,0.);
%    w1 = padarray( window(@tukeywin,(2520+wdSz(n)),50/(wdSz(n)+2520)), (szt(1) - (2520+wdSz(n)))/2.,0.);
      w1 = padarray( window(@tukeywin,(2520+wdSz(n)),0.), (szt(1) - (2520+wdSz(n)))/2.,0.);
  
    cropped = (truth.*w1).*w1';
    
    propagatedField = Frwrd* cropped;
    intensity = abs(propagatedField(rngt,rngt)).^2;
    szd = size(intensity);
    clear  propagatedField;
  %%  
    for m=1:numel(sr)
        S = LinOpDiag([sr(m)  szd(1)/sr(m)  sr(m) szd(2)/sr(m) nAngle],1./(sr(m).^2))*LinOpShape( [szd  ],[sr(m)  szd(1)/sr(m)  sr(m) szd(2)/sr(m) nAngle]);
        
        Sd = squeeze( sum( sum( S*trueIntensity,3),1));
        Sc = squeeze( sum( sum( S*intensity,3),1));
        szs = size(Sd);
        
        
        D = LinOpDiag([sr(m)  szt(1)/sr(m)  sr(m) szt(2)/sr(m) nAngle],1./(sr(m).^2))*LinOpShape( [szt  ],[sr(m)  szt(1)/sr(m)  sr(m) szt(2)/sr(m) nAngle]);
        dwn = squeeze( sum( sum( D*cropped,3),1));
        dwn = (LinOpPropagator(size(dwn),lambda, n0, z,dxy*sr(m),  Iangle, 1,1,'BLAS2'))*dwn;
        rngd = ceil((szt(1) - 2520)/sr(m)/2+1): ceil((szt(1) + 2520)/sr(m)/2);
        dwn = abs(dwn(rngd,rngd)).^2;
        %%
        dwn = dwn(:);
        Sc = Sc(:);
        
        parfor p=1:numel(SNR)
            tmpPSNR2 = 0.;
            tmpPSNR = 0.;
            tmpPSNR3 = 0.;
            for t=1:100
                tmpSd =  Sd + SNR(p).*random('Normal',zeros(szs),ones(szs));
                tmpPSNR2  = tmpPSNR2+ 10*log10(sum(abs(dwn-tmpSd(:)).^2)/(numel(Sd)-1));
                tmpPSNR  =  tmpPSNR +10*log10(sum(abs(Sc-tmpSd(:)).^2)/(numel(Sd)-1));
            end
            PSNR2(n,m,p) =  tmpPSNR2/100. ;
            PSNR(n,m,p) = tmpPSNR/100. ;
        end
        
        PSNR3(n,m) = 10*log10(sum(abs(dwn-Sc).^2)/(numel(Sd)-1));
        
    end
    clear cropped
    %dc = abs(yc(9501:10500,9501:10500)).^2;
    %MSE2 = [MSE2,20*log10(norm(data(:)-dc(:)))];
end