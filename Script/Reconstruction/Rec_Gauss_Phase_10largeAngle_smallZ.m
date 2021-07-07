useGPU(0);
%Simu_Phase_10largeAngle_smallZ;
%load('Simu_Phase_10largeAngle_SmallZ.mat');
% d1 = d(976:975+308,876:875+308,:);
%d1 = d(625:624+1024,525:524+1024,:);
%d1 = d(985:984+256,885:884+256,:);

load('Simu_Phase_10largeAngle_SmallZ_notilted.mat');

d1 = d(1026:1025+346,876:875+346,:);
%
% truth((9751-2060):(9750+3080+2060),(8751-2060):(8750+3080+2060))
%
%clear d
%% Physical parameters
k = n0*2*pi/lambda;
nAngle = size(d1,3);

%% Data generation
ill = 1; % Illumination
SNR = 60; 
sigma = 10^(-SNR./20);
eta = sigma./ill;
sizeData = size(d1);
data = ill.*d1;% + sigma.*random('Normal',zeros(sizeData),ones(sizeData));

%% Fov and badwidth
eta=0.035
outoffov = (sqrt( 2*z*k*pi*eta.^2*cos(theta)+ sin(theta).^2) + sin(abs(theta)) ) ./  (k *pi * eta.^2 * cos(theta).^2);

shift = z.*tan(theta);

pixsz = min(lambda./n0./ (outoffov./sqrt(outoffov.^2+z^2) + sin(abs(theta)))/2,[],2);
pixsz = min(pixsz) % square pixel
sr = ceil(dxy/pixsz);
sr=8
pixsz = dxy / sr;

fov = (max(outoffov+shift,[],2)+ max(outoffov-shift,[],2))' + sizeData(1:2)*dxy;
fovpix = ceil(ceil(fov/pixsz)/sr)*sr;


%sr = ceil(n0*dxy/(lambda/4));
%pixsz = dxy / sr;
%fovpix = ceil(ceil(fov/pixsz)/sr)*sr;
%lowfovpix = fovpix/sr;

% Bounding box
lbcorner = round(max(outoffov+shift,[],2)/dxy)';
leftext = lbcorner(1);
bottomext = lbcorner(2);

rtcorner = round(max(outoffov-shift,[],2)/dxy)';
rightext = rtcorner(1);
topext = rtcorner(2);

fovpix = (lbcorner +  rtcorner + sizeData(1:2))*sr%% Forward model
lowfovpix = fovpix/sr;

%%
% H = LinOpPropagator(fovpix,lambda, n0, z,pixsz,  theta, ill,1,'BLAS2','oversample');
%

%% Reshaping operator [szx  nAngle] -> [sr  szx(1)/sr  sr szx(2)/sr nAngle])
S = LinOpShape( [fovpix  nAngle],[sr  lowfovpix(1)  sr lowfovpix(2) nAngle]);
H = S* LinOpPropagatorH(fovpix,lambda, n0, z,pixsz,  theta,1./sr,1,'BLAS2');

%% Reshaping operator [szx  nAngle] -> [sr  szx(1)/sr  sr szx(2)/sr nAngle])
% 
% S = LinOpShape( [sizeData(1)*sr sizeData(2)*sr nAngle],[sr  sizeData(1)  sr sizeData(2) nAngle]);
% SL =  LinOpCrop([fovpix  nAngle],[ sr*lbcorner+1; sr*(-rtcorner)]);
% F =  LinOpPropagatorH(fovpix,lambda, n0, z,pixsz,  theta,1./sr,1,'BLAS2');
% H = S*SL*F;

%% Prior
D = LinOpHFGrad(fovpix);
rglCost = CostMixNorm21(D.sizeout,[3]);
C = LinOpSDFT(fovpix)^(-1);
constraint  = CostComplexCircle(fovpix,1.);

%% Likelihood
ldata =padarray(padarray(data/ill,[rightext,bottomext],1,'post'),[leftext,topext],1,'pre');
lwght =padarray(padarray( ones_(sizeData) ,[rightext,bottomext],'post'),[leftext,topext],'pre'); 
lklCost = CostIntensity([sr  lowfovpix(1)  sr lowfovpix(2) nAngle],ldata, lwght,[1 3],'Gaussian');
%lklCost = CostIntensity([sr  sizeData(1)  sr sizeData(2) nAngle],data, 1.,[1 3],'Gaussian');

clear lwght data d1

%% init
Hs = LinOpPropagator(lowfovpix,lambda, n0, z,dxy,  theta, 1,1,'BLAS2');
backprop = Hs'*ldata;
%xinit = constraint.applyProx(imresize(backprop,sr),1);
figure(3); imshow(angle(backprop),[1,2.5])
clear Hs  ldata


%% parameter
mu = 1e-11;
rho = 200*mu;
lklRho = rho
cRho  = 1.* rho*nAngle;
rglRho = rho/2*nAngle; 
Fn={lklCost,mu*rglCost,constraint};
Hn={H,D,C};
rho_n=[lklRho,rglRho,cRho];
ADMM=OptiADMM([],Fn,Hn,rho_n);
%ADMM.OutOp=OutputOptiSNR(1,im,10,[1 2]);
% STOP when the sum successives C = F*x + Fn{1}*Hn{1}*x is lower than 1e-4 or when the distance between two successive step is lower than 1e-5
%ADMM.CvOp=TestCvgCombine(TestCvgCostRelative(1e-4,[1 2]), 'StepRelative',1e-4);  
ADMM.ItUpOut=10;             % call OutputOpti update every ItUpOut iterations
ADMM.maxiter=1000;            % max number of iterations
ADMM.OutOp=OutputOptiCOMCI(false,10,true,[],C);
%%
%ADMM.run(fftn(constraint.applyProx(imresize(backprop,sr),1)));      % run the algorithm 
%ADMM.run(fftn(ones_( fovpix)));
ADMM.run(fftn(res));      % run the algorithm 
res = C*(ADMM.xopt);
