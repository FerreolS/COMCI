
PSpectCoherence = 0;

nbi=1;
load('TestSpectralCoherence_Gauss30.mat')
nbi = nb+1;
for nb = nbi:30
    nb
    switch randi(5,1)
        case 1
            SimuOrion;
        case 2
            SimuHelix;
        case 3 
            SimuBlueMarble;
        case 4
            SimuUSAFamp;
        case 5
            SimuUSAFph;
    end
           
   cshift = randi(10000,1,2)- 5000;
   truth= circshift(truth,cshift);
   TestSpectralCoherence;
   PSpectCoherence = PSpectCoherence + PSNR;
   save('TestSpectralCoherence_Gauss30.mat','nb','CoherenceLength','PSNR','IncidenceA','SNR','wdSz','lambda','k','z','n0','dxy','PSpectCoherence','dlambda');
end


PSpectCoherence = PSpectCoherence ./nb;