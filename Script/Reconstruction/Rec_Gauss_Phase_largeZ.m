run('Simu_Phase_largeZ.m');

useGPU(1);

nAngle = 1;
Iangle =[ 0 0]';
theta = Iangle./180*pi;
lambda= 532e-9; % wavelenght
n0 = 1.0;
k = 2.*pi/lambda*n0; % wavenumber
z = 265e-3;
 
dxy = 6.7e-6; % Camera  pixel size




%% Data generation
ill = 1; % Illumination
SNR = 35; 
sigma = 10^(-SNR./20);
eta = sigma./ill;
sizeData = size(d1);
data = ill.*d1 + sigma.*random('Normal',zeros(sizeData),ones(sizeData));

%% Fov and badwidth
outoffov = (sqrt( 2*z*k*pi*eta.^2*cos(theta)+ sin(theta).^2) + sin(abs(theta)) ) ./  (k *pi * eta.^2 * cos(theta).^2);
shift = z.*tan(theta);

pixsz = min(lambda./ (outoffov./sqrt(outoffov.^2+z^2) + sin(abs(theta)))/2,[],2);
pixsz = min(pixsz); % square pixel
sr = ceil(dxy/pixsz);
pixsz = dxy / sr;

fov = ((max(outoffov,[],2) + max( abs(shift),[],2))'.*2 + sizeData(1:2)*dxy);
fovpix = ceil(ceil(fov/pixsz)/sr)*sr;

%% Bounding box
lbcorner = round(max(outoffov+shift,[],2)/dxy)';
leftext = lbcorner(1);
bottomext = lbcorner(2);

rtcorner = round(max(outoffov-shift,[],2)/dxy)';
rightext = rtcorner(1);
topext = rtcorner(2);

fovpix = lbcorner +  rtcorner + sizeData(1:2)*sr;
%% Forward model

lowfovpix = fovpix/sr;
Hs = LinOpPropagator(lowfovpix,lambda, n0, z,dxy,  theta, ill,1,'BLAS2');

%% Reshaping operator [szx  nAngle] -> [sr  szx(1)/sr  sr szx(2)/sr nAngle])
S = LinOpShape( [fovpix  nAngle],[sr  lowfovpix(1)  sr lowfovpix(2) nAngle]);
H = S* LinOpPropagatorH(fovpix,lambda, n0, z,pixsz,  theta, sqrt(ill)./sr,1,'BLAS2');

%% Likelihood
ldata =padarray(padarray(data,[rightext,bottomext],0.,'post'),[leftext,topext],0.,'pre');
lwght =padarray(padarray( ones_(sizeData) ,[rightext,bottomext],'post'),[leftext,topext],'pre'); 
gpuCpuConverter(ldata);
gpuCpuConverter(lwght);

lklCost = CostIntensity([sr  lowfovpix(1)  sr lowfovpix(2) nAngle],ldata, lwght,[1 3],'Gaussian');


%% Prior
D = LinOpHFGrad(fovpix);
rglCost = CostMixNorm21(D.sizeout,[3]);
C = LinOpSDFT(fovpix)^(-1);
constraint  = CostComplexCircle(fovpix,1.);

%% init
backprop = Hs'*ldata;
xinit = constraint.applyProx(imresize(abs(backprop),sr),1);
xinit = ones_(size(xinit));

%% parameter
mu = 10.^-1;
rho = mu.*10.^(2);
lklRho = rho/nAngle;
cRho  = 1.* rho;
rglRho = rho;

%% -- ADMM LS + TV + NonNeg
Fn={lklCost,mu*rglCost,constraint};
Hn={H,D,C};
rho_n=[lklRho,rglRho,cRho];
ADMM=OptiADMM([],Fn,Hn,rho_n);
%ADMM.OutOp=OutputOptiSNR(1,im,10,[1 2]);
% STOP when the sum successives C = F*x + Fn{1}*Hn{1}*x is lower than 1e-4 or when the distance between two successive step is lower than 1e-5
%ADMM.CvOp=TestCvgCombine(TestCvgCostRelative(1e-4,[1 2]), 'StepRelative',1e-4);  
ADMM.ItUpOut=10;             % call OutputOpti update every ItUpOut iterations
ADMM.maxiter=1000;            % max number of iterations
ADMM.run(fft2(xinit));      % run the algorithm 
res = C*(ADMM.xopt);
