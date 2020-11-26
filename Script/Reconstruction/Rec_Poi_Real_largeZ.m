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

sizeData = size(d1);
nph=3;
Shot = 100;
ill=10.^(nph);
dph = ill.* gather(d1);
data = random('Poisson',dph + Shot);

sigma = std( data(:) - (dph(:) +Shot(:)));
clear dph
eta = gather(sigma./ill);
%eta = sqrt(Shot)/ill;

%%
%I1p = 10.^(nph/2.+2.)./sum(I1(:)) .* I1;
%I1n = random('Poisson',I1p + Shot);

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

fovpix = (lbcorner +  rtcorner + sizeData(1:2))*sr;
%% Forward model
lowfovpix = fovpix/sr;
Hs = LinOpPropagator(lowfovpix,lambda, n0, z,dxy,  theta, sqrt(ill),1,'BLAS2');

%% Reshaping operator [szx  nAngle] -> [sr  szx(1)/sr  sr szx(2)/sr nAngle])
S = LinOpShape( [fovpix*2  nAngle],[sr*2  lowfovpix(1)  sr*2 lowfovpix(2) nAngle]);
H = S* LinOpPropagatorH(fovpix,lambda, n0, z,pixsz,  theta, sqrt(ill)./sr,1,'BLAS2','oversample');

%% Likelihood
ldata =padarray(padarray(data,double([rightext,bottomext]),0.,'post'),double([leftext,topext]),0.,'pre');
lwght =padarray(padarray(Shot.*ones_(sizeData) ,[rightext,bottomext],-1,'post'),[leftext,topext],-1,'pre'); 

gpuCpuConverter(ldata);
gpuCpuConverter(lwght);

lklCost = CostIntensity([sr*2  lowfovpix(1)  sr*2 lowfovpix(2) nAngle],ldata, lwght,[1 3],'Poisson');

%lwght =padarray(padarray(1./max(Shot,data) ,[rightext,bottomext],0,'post'),[leftext,topext],0,'pre'); 
%gpuCpuConverter(lwght);
%lklCost = CostIntensity([sr  lowfovpix(1)  sr lowfovpix(2) nAngle],ldata, lwght,[1 3],'Gaussian');


%% Prior
D = LinOpHFGrad(fovpix);
rglCost = CostMixNorm21(D.sizeout,[3]);
C = LinOpSDFT(fovpix)^(-1);
constraint  = CostReals(fovpix,0.,1.);


%% init
backprop = Hs'*ldata;
xinit = constraint.applyProx(imresize(abs(backprop),sr),1);
xinit = ones_(size(xinit));

%% parameter
mu = 1e-1;
rho = mu.*10.^(2);
lklRho = rho/sqrt(ill)* prod(2*fovpix)/numel(data);
cRho  = 1.* rho;
rglRho = rho/2.;

% -- ADMM LS + TV + NonNeg
Fn={lklCost,mu*rglCost,constraint};
Hn={H,D,C};
rho_n=[lklRho,rglRho,cRho];
ADMM=OptiADMM([],Fn,Hn,rho_n);
%ADMM.OutOp=OutputOptiSNR(1,im,10,[1 2]);
% STOP when the sum successives C = F*x + Fn{1}*Hn{1}*x is lower than 1e-4 or when the distance between two successive step is lower than 1e-5
%ADMM.CvOp=TestCvgCombine(TestCvgCostRelative(1e-4,[1 2]), 'StepRelative',1e-4);  
ADMM.ItUpOut=10;             % call OutputOpti update every ItUpOut iterations
ADMM.maxiter=100;            % max number of iterations
ADMM.run(fft2(xinit));      % run the algorithm 
res = C*(ADMM.xopt);
figure; imshow(abs(res),[])
