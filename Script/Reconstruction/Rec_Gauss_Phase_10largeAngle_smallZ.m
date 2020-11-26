
load('Angle_AS_Simu.mat');
d1 = d(1181:1180+512,1181:1180+512,[1:7,9,10]);
% d1 = d(1310:1309+256,1310:1309+256,[1:7,9,10]);

 clear d
%% Physical parameters
dxy = 1.12e-6;   % pixelsize
lambda= 700e-9; % wavelenght
z = 100e-6 ;  % depth %
n0 = 1.0;       % refractive index of the propagation medium
k = 2.*pi/lambda*n0; % wavenumber
IangleT =[ [  -50    -25     0     25     50     0     0         0 0];[   0     0     0     0     0     -50     -25       25 50] ];
theta = IangleT./180*pi;
nAngle = size(d1,3);

%% Data generation
ill = 1; % Illumination
SNR = 25; 
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

fov = (max(outoffov+shift,[],2)+ max(outoffov-shift,[],2))' + sizeData(1:2)*dxy;
fovpix = ceil(ceil(fov/pixsz)/sr)*sr;


%sr = ceil(dxy/(lambda/4));
%pixsz = dxy / sr;
%fovpix = ceil(ceil(fov/pixsz)/sr)*sr;
%lowfovpix = fovpix/sr;

%% Bounding box
lbcorner = round(max(outoffov+shift,[],2)/dxy)';
leftext = lbcorner(1);
bottomext = lbcorner(2);

rtcorner = round(max(outoffov-shift,[],2)/dxy)';
rightext = rtcorner(1);
topext = rtcorner(2);

fovpix = (lbcorner +  rtcorner + sizeData(1:2))*sr;%% Forward model
lowfovpix = fovpix/sr;

%%
% H = LinOpPropagator(fovpix,lambda, n0, z,pixsz,  theta, ill,1,'BLAS2','oversample');
%

%% Reshaping operator [szx  nAngle] -> [sr  szx(1)/sr  sr szx(2)/sr nAngle])
S = LinOpShape( [fovpix  nAngle],[sr  lowfovpix(1)  sr lowfovpix(2) nAngle]);
H = S* LinOpPropagatorH(fovpix,lambda, n0, z,pixsz,  theta, sqrt(ill)./sr,1,'BLAS2');

%% Prior
D = LinOpHFGrad(fovpix);
rglCost = CostMixNorm21(D.sizeout,[3]);
C = LinOpSDFT(fovpix)^(-1);
constraint  = CostComplexCircle(fovpix,1.);

%% Likelihood
ldata =padarray(padarray(data,[rightext,bottomext],0.,'post'),[leftext,topext],0.,'pre');
lwght =padarray(padarray( sigma.^(-2).*ones_(sizeData) ,[rightext,bottomext],'post'),[leftext,topext],'pre'); 
lklCost = CostIntensity([sr  lowfovpix(1)  sr lowfovpix(2) nAngle],ldata, lwght,[1 3],'Gaussian');



%% init
% Hs = LinOpPropagator(lowfovpix,lambda, n0, z,dxy,  theta, ill,1,'BLAS2');
% backprop = Hs'*ldata;
% xinit = constraint.applyProx(imresize(backprop,sr),1);
% clear Hs


%% parameter
mu = 1e-4;
rho = mu.*10.^(2);
lklRho = rho/nAngle/(ill)*sigma.^2;
cRho  = 1.* rho;
rglRho = rho/2; 
Fn={lklCost,mu*rglCost,constraint};
Hn={H,D,C};
rho_n=[lklRho,rglRho,cRho];
ADMM=OptiADMM([],Fn,Hn,rho_n);
%ADMM.OutOp=OutputOptiSNR(1,im,10,[1 2]);
% STOP when the sum successives C = F*x + Fn{1}*Hn{1}*x is lower than 1e-4 or when the distance between two successive step is lower than 1e-5
%ADMM.CvOp=TestCvgCombine(TestCvgCostRelative(1e-4,[1 2]), 'StepRelative',1e-4);  
ADMM.ItUpOut=2;             % call OutputOpti update every ItUpOut iterations
ADMM.maxiter=100;            % max number of iterations
ADMM.run(fftn(ones_(fovpix)));      % run the algorithm 
res = C*(ADMM.xopt);
