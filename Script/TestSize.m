szt = size(truth);

nAngle = 1;
Iangle = zeros(2,nAngle);
lambda= 530e-9; % wavelenght
n0 = 1.0;
z = 250*1e-6;
%dxy = 5.2e-6; % Camera  pixel size
dxy= 1.04e-6 ;
sr=16; % dxy / dxyin


Frwrd = LinOpPropagator(szt,lambda, n0, z,dxy/sr,  Iangle, 1,1,'AS');
y = Frwrd* truth;

data = abs(y(9501:10500,9501:10500)).^2;
MSE = [];
for n=1:19
    newSz = 500*n;
    cropped = truth(1+newSz:szt(1)-newSz,1+newSz:szt(1)-newSz);
    szc = size(cropped);
    Frwrd = LinOpPropagator(szc,lambda, n0, z,dxy/sr,  Iangle, 1,1,'AS');
    yc = Frwrd* cropped;
    dc = abs(yc(9501-newSz:10500-newSz,9501-newSz:10500-newSz)).^2;
    MSE = [MSE,20*log10(norm(data(:)-dc(:)))];
end