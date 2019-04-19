szt = size(truth);

nAngle = 1;
Iangle = zeros(2,nAngle);
lambda= 530e-9; % wavelenght
n0 = 1.0;
z = 250*1e-6;
k = n0*2.*pi/lambda;
%dxy = 5.2e-6; % Camera  pixel size
dxy= 1.04e-6 / 10;


sr=[1:7 9 10 12 14 15 18 20];  %downsampling factor
SNR = [0 0.001 0.002 0.005 0.0075 0.01 0.015 0.02 0.05 0.075 0.1 0.2 0.5 1] ;
wdSz = [0 0.01 0.02 0.05 0.1 0.15 0.2 0.3 0.4 0.5 0.75 1. 1.25 1.5  1.75 2 3 4 5 6 7]*1000;
IncidenceA = [-60 -45 -30 -15 0 15 30 45 60];

PSNR = zeros(numel(wdSz),numel(IncidenceA),numel(SNR));
m=1;

%%
for a=1:numel(IncidenceA)
    nominalAngle = [IncidenceA(a)/180*pi; 0];
    w = padarray( window(@tukeywin,szt(1)/2,0.), szt(1)/4,0.);
    truth = (truth.*w).*w';
    
    Frwrd = LinOpPropagator(szt,lambda, n0, z,dxy,  nominalAngle, 1,1,'BLAS2');
    truePropagatedField = Frwrd* (truth);
    
    rngt = ((szt(1) - 2520)/2+1): (szt(1) + 2520)/2;
    trueIntensity = abs(truePropagatedField(rngt,rngt)).^2;
    %%
    clear truePropagatedField;
    for n=1:numel(wdSz)
        %%  w1 = padarray( window(@tukeywin,1000*n,1000*n/400), (20000 - 1000*n)/2.,0.);
        %    w1 = padarray( window(@tukeywin,(2520+wdSz(n)),50/(wdSz(n)+2520)), (szt(1) - (2520+wdSz(n)))/2.,0.);
        %     w1 = padarray( window(@tukeywin,(2520+wdSz(n)),0.), (szt(1) - (2520+wdSz(n)))/2.,0.);
        %    cropped = (truth.*w1).*w1';
        cropped = truth;
        cropped(((szt(1) + 2520)/2 +wdSz(n)):end,:) =0;
        
        propagatedField = Frwrd* cropped;
        intensity = abs(propagatedField(rngt,rngt)).^2;
        szd = size(intensity);
        clear  propagatedField;
     
       
        Sc = intensity(:);
        Sd = trueIntensity(:);
        szs = size(Sd);
        
        parfor p=1:numel(SNR)
            tmpPSNR = 0.;
            for t=1:100
                tmpSd =  Sd + SNR(p).*random('Normal',zeros(szs),ones(szs));
                tmpPSNR  =  tmpPSNR +10*log10(sum(abs(Sc-tmpSd(:)).^2)/(numel(Sd)-1));
            end
            PSNR(n,a,p) = tmpPSNR/100. ;
        end
                
    end
    clear cropped
end