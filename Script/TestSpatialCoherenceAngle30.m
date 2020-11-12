
PSpatCoherence = 0;

nbi=1;
load('TestSpatialCoherence_Disk_Angle30.mat')
nbi = nb;
%
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
   TestSpatialCoherence;
   PSpatCoherence = PSpatCoherence + PSNR;
   save('TestSpatialCoherence_Disk_Angle30.mat','HalfApexAngle','nb','PSNR','IncidenceA','SNR','wdSz','lambda','k','z','n0','dxy','PSpatCoherence');
end


PSpatCoherence = PSpatCoherence ./nb;