szt = size(truth);

w = padarray( window(@tukeywin,szt(1)/2,0.), szt(1)/4,0.);
truth = (truth.*w).*w';
%%
lambda= 530e-9; % wavelenght
n0 = 1.0;
z = 250*1e-6;
k = n0*2.*pi/lambda;
%dxy = 5.2e-6; % Camera  pixel size
dxy= 1.04e-6 / 10;
c = 299792458; % speed of ligth
dlambda = 10.e-9; % line fwhm
dsig = dlambda/(2*sqrt(2*log(2)));
    N = 100;
    
%%
sr=[1:7 9 10 12 14 15 18 20];  %downsampling factor
SNR = [0 0.001 0.002 0.005 0.0075 0.01 0.015 0.02 0.05 0.075 0.1 0.2 0.5 1] ;
wdSz = [0  0.05 0.1 0.3 0.5 0.7 0.9 1. 1.1 1.3 1.5 1.7 1.9 2.1 2.25 2.5 3 3.5 4]*1000;
IncidenceA = [ -45 -30 -15 0 15 30 45];
PSNR = zeros(numel(wdSz),numel(IncidenceA),numel(SNR));

%%
for a=1:numel(IncidenceA)
    %%
    nAngle = 1;
    nominalAngle = [IncidenceA(a)/180*pi; 0];
    
    
    trueIntensity=0;
    centerx = -round(z/dxy* tan(nominalAngle(1)));
    rngtx = (centerx+(szt(1) - 2520)/2+1): (centerx+(szt(1) + 2520)/2);
    rngty = ((szt(1) - 2520)/2+1): (szt(1) + 2520)/2;
    Frwrd = LinOpPropagator(szt,lambda, n0, z,dxy, nominalAngle, 1,1,'BLAS2');
   % truePropagatedField = Frwrd* truth;
   % coherentIntensity=abs(truePropagatedField(rngtx,rngty)).^2;
    clear  truePropagatedField;
    
    %%
    for tt=1:N
        lmbd = random( 'norm',0.,dsig) + lambda;
        Frwrd.lambda =   lmbd;
        truePropagatedField = Frwrd* truth;
        trueIntensity =trueIntensity+ 1./N*abs(truePropagatedField(rngtx,rngty)).^2;
        clear  truePropagatedField;
    end
    %%
    Frwrd.lambda =   lambda;
    
    for n=1:numel(wdSz)
        cropped = truth;
        cropped(((szt(1) + 2520)/2 +wdSz(n)):end,:) =0;
        
        propagatedField = Frwrd* cropped;
        intensity = abs(propagatedField(rngtx,rngty)).^2;
        szd = size(intensity);
        clear  propagatedField;
        
        Sc = intensity(:);
        Sd = trueIntensity(:);
        szs = size(Sd);
            
        parfor p=1:numel(SNR)
            tmpPSNR = 0.;
            for t=1:10
                tmpSd =  Sd + SNR(p).*random('Normal',zeros(szs),ones(szs));
                tmpPSNR  =  tmpPSNR +10*log10(sum(abs(Sc-tmpSd(:)).^2)/(numel(Sd)-1));
            end
            PSNR(n,a,p) = tmpPSNR/10. ;
        end
        
        clear cropped
        
    end
end