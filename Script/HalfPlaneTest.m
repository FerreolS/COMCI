
lambda= 530e-9; % wavelenght
n0 = 1.0;
z = 250*1e-6;
k = n0*2.*pi/lambda;
%dxy = 5.2e-6; % Camera  pixel size
dxy= 1.04e-6 / 16;

c = sqrt(k/pi/z).*((1:2520) + 500)*dxy;
truth = repmat(1/2 * ( (1/2 + fresnelc(c)).^2 +  (1/2 + fresnels(c)).^2),2520,1);
szt = size(truth);

trueIntensity = truth+10.*random('Normal',zeros(szt),ones(szt));

SNR = [0 0.001 0.002 0.005 0.0075 0.01 0.015 0.02 0.05 0.075 0.1 0.2 0.5 1] ;

for n=1:numel(SNR)
    trueIntensity = truth+ SNR(n).*random('Normal',zeros(szt),ones(szt));
