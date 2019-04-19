P1 = 0;
P2 = 0;
P3 = 0;

load('TestSize2_x100.mat');
nbi = nb+1;
%nbi=1;
for nb = nbi:80
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
   TestSize2;
   P1 = P1 + PSNR;
   P2 = P2 + PSNR2;
   P3 = P3 + PSNR3;
   save('TestSize2_x100.mat','nb','P1','P2','P3', 'PSNR','PSNR2','PSNR3','sr','SNR','wdSz','Iangle','lambda','k','z','n0','dxy');
end


P1 = P1 ./nb;
P2 = P2 ./nb;
P3 = P3 ./nb;