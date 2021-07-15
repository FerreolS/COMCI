
%%
clear all
useGPU(0)
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



H = LinOpPropagatorH(fovpix,lambda, n0, z,pixsz,  theta, sqrt(ill)./sr,1,'Fresnel');


%% Likelihood

ldata =padarray(padarray(data,[rightext,bottomext],0.,'post'),[leftext,topext],0.,'pre');
lwght =padarray(padarray( ones_(sizeData) ,[rightext,bottomext],'post'),[leftext,topext],'pre'); 


%% parameter
mu = 10.^-2;
xinit = ones_(fovpix);


%% -- VMLMB LS + hyperbolicTV + non negativity
% H = LinOpPropagatorH(fovpix,lambda, n0, z,pixsz,  theta, sqrt(ill)./sr,1,'Fresnel');
H = LinOpPropagator(fovpix,lambda, n0, z,pixsz,  theta, sqrt(ill)./sr,1,'Fresnel');

M = OpEWSquaredMagnitude(fovpix);
W = LinOpDiag(fovpix,lwght);
LS=CostL2(fovpix,ldata,W)*M;                 % Least-Squares data term

C = LinOpCpx(fovpix); % utility operator to go from complex to reals pair
F = LinOpSDFT(fovpix)^(-1);
D = LinOpGrad(fovpix);
hyperB = CostHyperBolic(D.sizeout,   1e-3,  3)*D;
%hyperB = CostL2(D.sizeout)*D;
L = (LS*H+ mu*hyperB)*C'; 
L.memoizeOpts.apply=true;
VMLMBPos=OptiVMLMB(L,0,[]);  
VMLMBPos.OutOp=OutputOptiSNR(1,C*t,1,[1 2]);
%VMLMBPos.CvOp=TestCvgCombine('CostRelative',1e-6, 'StepRelative',1e-6);
VMLMBPos.ItUpOut=1; 
VMLMBPos.maxiter=1000;                             % max number of iterations
VMLMBPos.m=5;                                     % number of memorized step in hessian approximation

VMLMBPos.run(C*xinit);                                  % run the algorithm 
resPos = C'*(VMLMBPos.xopt);
%save('ExtLinearO30dB_VMLM_Pos.mat','VMLMBPos');

%% -- VMLMB LS + hyperbolicTV
VMLMB=OptiVMLMB(L,[],[]);  
VMLMB.OutOp=OutputOptiSNR(1,C*t,1,[1 2]);
%VMLMB.CvOp=TestCvgCombine('CostRelative',1e-6, 'StepRelative',1e-6);
VMLMB.ItUpOut=1; 
VMLMB.maxiter=1000;                             % max number of iterations
VMLMB.m=5;                                     % number of memorized step in hessian approximation

%%
VMLMB.run(C*xinit);                                  % run the algorithm 
res = C'*(VMLMB.xopt);
%save('ExtLinearO30dB_VMLM.mat','VMLMB');
