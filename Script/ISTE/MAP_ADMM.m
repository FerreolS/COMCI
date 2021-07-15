
%%
clear all
useGPU(1)
SimulateModelFieldExtension
sizeData = size(IntensityModel);
%% SNR = 30
ill = 1; % Illumination
SNR = 30;
sigma = 10^(-SNR./20).*ill;
data = max(0.,IntensityModel.*ill+ sigma.*(random('Normal',zeros(sizeData),ones(sizeData)))) ;
%t = truth(1537:1536+1024,1537:1536+1024,:);
t = truth(973:972+2152,973:972+2152,:);


%%

eta = sigma./ill;
outoffov = (sqrt( 2*z*k*pi*eta.^2*cos(theta)+ sin(theta).^2) + sin(abs(theta)) ) ./  (k *pi * eta.^2 * cos(theta).^2);

sr =1;
pixsz = dxy / sr;
fov = ((max(outoffov,[],2))'.*2 + sizeData(1:2)*dxy);
fovpix = ceil(ceil(fov/pixsz)/sr)*sr;

disp(['the fov extrapolation is ',num2str(outoffov(1)./dxy), ' pixels in each direction' ])
disp(['leading to fov of width ',num2str(2*outoffov(1)./dxy+1024), ' pixels ( ', num2str((2*outoffov(1)./dxy+1024)*dxy*1e3), ' mm)'])

%% Bounding box
lbcorner = round(max(outoffov,[],2)/dxy)';
leftext = lbcorner(1);
bottomext = lbcorner(2);

rtcorner = round(max(outoffov,[],2)/dxy)';
rightext = rtcorner(1);
topext = rtcorner(2);

fovpix = (lbcorner +  rtcorner + sizeData(1:2))*sr;
%% Forward model

lowfovpix = fovpix/sr;
Hs = LinOpPropagator(lowfovpix,lambda, n0, z,dxy,  theta, ill,1,'BLAS2');


%% Reshaping operator [szx  nAngle] -> [sr  szx(1)/sr  sr szx(2)/sr nAngle])
S = LinOpShape( [fovpix*2  1],[sr*2  lowfovpix(1)  sr*2 lowfovpix(2) 1]);
H = S* LinOpPropagatorH(fovpix,lambda, n0, z,pixsz,  theta, sqrt(ill)./sr,1,'Fresnel','oversample');


%% Likelihood

ldata =padarray(padarray(data,[rightext,bottomext],0.,'post'),[leftext,topext],0.,'pre');
lwght =padarray(padarray( ones_(sizeData) ,[rightext,bottomext],'post'),[leftext,topext],'pre'); 
gpuCpuConverter(ldata);
gpuCpuConverter(lwght);

lklCost = CostIntensity([sr*2  lowfovpix(1)  sr*2 lowfovpix(2) nAngle],ldata, lwght,[1 3],'Gaussian');

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
mu = 10.^-2;
rho = mu.*10.^(2);
lklRho = rho/nAngle;
cRho  = 1.* rho;
rglRho = rho;

%% -- ADMM  TV 
Fn={lklCost,mu*rglCost,constraint};
Hn={H,D,C};
rho_n=[lklRho,rglRho,cRho];
ADMM=OptiADMM([],Fn,Hn,rho_n);
ADMM.OutOp=OutputOptiSNR(1,C^(-1)*t,1,[1 2]);
% STOP when the sum successives C = F*x + Fn{1}*Hn{1}*x is lower than 1e-4 or when the distance between two successive step is lower than 1e-5
%ADMM.CvOp=TestCvgCombine(TestCvgCostRelative(1e-4,[1 2]), 'StepRelative',1e-4);  
ADMM.ItUpOut=1;             % call OutputOpti update every ItUpOut iterations
ADMM.maxiter=1000;            % max number of iterations
%%
ADMM.run(C^(-1)*xinit);      % run the algorithm 
res = C*(ADMM.xopt);

