P1 = 0;

nbi=1;
%load('TestSize_x100.mat')
%nbi = nb+1;
for nb = nbi:100
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
   TestSize;
   P1 = P1 + PSNR;
   save('TestSize_x100.mat','nb','P1', 'PSNR','wdSz','SNR','IncidenceA','lambda','k','z','n0','dxy');
send


P1 = P1 ./nb;