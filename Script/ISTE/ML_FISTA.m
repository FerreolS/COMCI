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

%% Likelihood

ldata =padarray(padarray(data,[rightext,bottomext],0.,'post'),[leftext,topext],0.,'pre');
lwght =padarray(padarray( ones_(sizeData) ,[rightext,bottomext],'post'),[leftext,topext],'pre'); 
gpuCpuConverter(ldata);
gpuCpuConverter(lwght);



%% parameter
mu = 10.^-2;
xinit = ones_(fovpix);




%% Operators
H = LinOpPropagator(fovpix,lambda, n0, z,pixsz,  theta, sqrt(ill)./sr,1,'Fresnel');
M = OpEWSquaredMagnitude(fovpix);
W = LinOpDiag(fovpix,lwght);

LS=CostL2(fovpix,ldata,W)*M;                 % Least-Squares data term

lklCost = CostIntensity([  fovpix(1) fovpix(2) ],ldata, lwght,'Gaussian');

%% -- FISTA LS + Phase only
Pphaseonly = CostComplexCircle(fovpix,1);

F=LS*H;
F.doPrecomputation=1;

FBS=OptiFBS(F,Pphaseonly);
FBS.gam=1.;
FBS.updateGam='backtracking';
FBS.OutOp=OutputOptiSNR(1,t,1,[1 2]);
%FBS.CvOp=TestCvgCombine(TestCvgCostRelative(1e-4), 'StepRelative',1e-4);  
FBS.ItUpOut=1;          % call OutputOpti update every ItUpOut iterations
FBS.fista=true;         % activate fista
FBS.maxiter= ;        % max number of iterations
FBS.run(xinit);% run the algorithm (Note that gam is fixed automatically to 1/F.lip here since F.lip is defined and since we do not have setted gam) 


%% -- Douglas-Rachford LS + Phase only

DR=OptiDouglasRachford(Pphaseonly,lklCost,H,1,1.);
DR.OutOp=OutputOptiSNR(1,t,1,[1 2]);
%DR.CvOp=TestCvgCombine(TestCvgCostRelative(1e-8), 'StepRelative',1e-8);  
DR.ItUpOut=1;          % call OutputOpti update every ItUpOut iterations
DR.maxiter=26;        % max number of iterations
DR.run(xinit);

%% --  ER LS + Phase Only
o = xinit;
maxiter = 100;
E = zeros_(maxiter,1);
S = zeros_(maxiter,1);
normXtrue=norm(t(:));
for n=1:maxiter
    dd = H*o;
    d = lklCost.applyProx(dd,1);
    oo = H'*d;
    o = Pphaseonly.applyProx(oo,1);
    E(n) = 20.*log10(norm(o(:)-oo(:)) +norm(d(:)-dd(:)));
    S(n) = 20.*log10(normXtrue/norm(t(:)-o(:)));
    disp(['iter : ', num2str(n), ' error :', num2str(E(n)), ' SNR :',num2str(S(n))])
end


%% --  ER strict + Phase Only
lklstrictCost = CostIntensity([  fovpix(1) fovpix(2) ],ldata, lwght,'Strict');

o = xinit;
maxiter = 120;
E = zeros_(maxiter,1);
S = zeros_(maxiter,1);
normXtrue=norm(t(:));
for n=1:maxiter
    dd = H*o;
    d = lklstrictCost.applyProx(dd,1);
    oo = H'*d;
    o = Pphaseonly.applyProx(oo,1);
    E(n) = 20.*log10(norm(o(:)-oo(:)) +norm(d(:)-dd(:)));
    S(n) = 20.*log10(normXtrue/norm(t(:)-o(:)));
    disp(['iter : ', num2str(n), ' error :', num2str(E(n)), ' SNR :',num2str(S(n))])
end

