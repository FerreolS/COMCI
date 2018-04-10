Folder = '/Users/ferreol/Google Drive/Pour Ferreol sans prisme/';
FolderShift = '/Users/ferreol/Google Drive/Pour Ferreol sans prisme/shifts/';
nAngle = 9;
theta = [-9.5 -9.5 -9.5 0 0 0 9.5 9.5 9.5];
phi = [-9.5 0  9.5 -9.5 0 9.5 -9.5 0  9.5];

Iangle = zeros(2,nAngle);
% 
% data = zeros(1024,1280,nAngle);
%  bg = zeros(1024,1280,nAngle);
%  rangex =1:1024;
%  rangey=1:1280;
x0 = 261;
y0 = 601;
Nx  = 340;
Ny = 340;
% 
% x0 = 1;
% y0 = 201;
% Nx  = 940;
% Ny = 940;
% 
% x0 = 1;
% y0 = 201;
% Nx  = 1024;
% Ny = 1024;

data = zeros(Nx,Ny,nAngle);
 bg = zeros(Nx,Ny,nAngle);
 rangex=x0:(x0+Nx-1);
 rangey =y0:(y0+Ny-1);
 
H = cell(nAngle,1);
flux =  zeros(1,nAngle);

dxy = 5.2e-6;
lambda= 681e-9; % wavelenght
z = 1.0256e-3;
n0 = 1.0;
eta = 0.01;
%%
sigma = 0.1;
ill=1;
Imat = zeros(2,nAngle);
for nA = 1:nAngle
    nt = theta(nA);
    np = phi(nA);
    Imat(:,nA)    = xlsread([FolderShift,'holo_shifts_',num2str(nt),'_',num2str(np),'.xls']);
    Iangle(1,nA) = -atan(Imat(1,nA)*dxy/z);
    Iangle(2,nA) = -atan(Imat(2,nA)*dxy/z); 
    file_path = [Folder,'holo_',num2str(nt),'_',num2str(np),'_bg.tif'];
    img = double(importdata(file_path));
   img =  img(:,:,1);
   img = img./ mean( img(:));
    tmpb= img(rangex,rangey);
    bg(:,:,nA) = tmpb;
    file_path = [Folder,'holo_',num2str(nt),'_',num2str(np),'.tif'];
    img = double(importdata(file_path));
   % img = img ./imgaussfilt(img,5);
   img =  img(:,:,1);
   img = img./ mean( img(:));
    img = img(rangex,rangey);
    
  % img = (img - min(img(:)))./(  max(img(:)) - min(img(:)));
    data(:,:,nA) = img ./tmpb  ;
  %  data(:,:,nA) = img ;
    flux(nA) =  mean(data(:,nA));%.* prod(cos( iangle),2);
end
%%

Iangle =[[ -0.080942278439585   0.059740376170844   0.201673699886618  -0.143317087312777   0.000118893052648   0.143921975794711  -0.201891376414567  -0.058494350395829   0.083751053815474]
   [0.193698187263600   0.133329522531104   0.073670203678171   0.059027710957567  -0.000097599727386  -0.061617864122155  -0.080153026577227  -0.142386941441917  -0.202347664250446]];

  
z =     0.001189416698366;

%%

sizeData = size(data(:,:,1));
Nx = sizeData(1);
Ny = sizeData(2);


lambdaT = lambda*ones([1, nAngle]);
 zT = z*ones([1, nAngle]);

 
 
% Bounds on angle
maxangle = max(Iangle,[],2);
minangle = min(Iangle,[],2);

%Bounds on sampled frequencies

asrx = 2.*max(abs(sin(Iangle(1,:)) ./  lambdaT))*  dxy  + 1;
asry = 2.*max(abs(sin(Iangle(2,:)) ./  lambdaT))*  dxy  + 1;

srx = ceil(asrx) ;
sry = ceil(asry) ;
sr = (max(srx,sry))
sr=4;

% Bounds on reconstructed object
tmp = ceil(abs( z*tan(minangle)/ dxy))+ ceil( (dxy/sr * max(lambdaT(:) .*zT(:)) /(8*pi*pi*eta^2)).^(1/3)/dxy);
leftext = tmp(1);
bottomext = tmp(2);
xlowcorner = tmp;
tmp = (ceil(abs( z*tan(maxangle))/ dxy))+ ceil( (dxy/sr * max(lambdaT(:) .*zT(:)) /(8*pi*pi*eta^2)).^(1/3)/dxy);
rightext = tmp(1);
topext = tmp(2);
xhighcorner = tmp;


% Reconstructed fov at camera resolution
fovlow = [( Nx+ leftext +  rightext),( Ny+ topext +  bottomext)];
bestfovlow = fft_best_dim(fovlow );
rightext = rightext + ceil((bestfovlow(1) - fovlow(1))/2);
leftext = leftext + floor((bestfovlow(1) - fovlow(1))/2);

topext = topext +ceil((bestfovlow(2) - fovlow(2))/2);
bottomext = bottomext +floor((bestfovlow(2) - fovlow(2))/2);
fovlow = bestfovlow;


szx =  sr.*fovlow ; % Size of the reconstructed object
szy = [ sr.* Nx,  sr.* Ny,  nAngle];


xleftext =  leftext* sr;
xbottomext =  bottomext* sr;
xrightext =  rightext* sr;
xtopext =  topext* sr;

%%


ldata =padarray(padarray(data,[rightext,bottomext],'replicate','post'),[leftext,topext],'replicate','pre');
lbg =padarray(padarray(bg,[rightext,bottomext],'replicate','post'),[leftext,topext],'replicate','pre');
lwght =padarray(padarray(ones(size(data))./max(sigma.^2,0.001),[rightext,bottomext],'post'),[leftext,topext],'pre');

lklCost = CostIntensity([sr  szx(1)/sr  sr szx(2)/sr nAngle],ldata , lwght,[1 3] ,'Gaussian');

%%
S = LinOpShape( [szx  nAngle],[sr  szx(1)/sr  sr szx(2)/sr nAngle]);

F = S*LinOpPropagatorH(szx,lambdaT, n0, zT,dxy/sr,  Iangle , sqrt(ill)./sr,'AS');

%%
 D = LinOpHFGrad(szx);
rglCost = CostMixNorm21(D.sizeout,[3,4]);


C =( LinOpDFT(szx))^(-1);
cCost  = CostComplexCircle(szx,1.);

%%
xinit =  ones(szx) ;

 
clear data bg  lbg lwght img
 %%   
mu  =5.;%2.5;

rho =    mu*100.;
 lklRho= rho/nAngle/ill;
rglRho = rho;
cRho  = 1.* rho;
%%

lklCost.doPrecomputation=1;
lklCost.memoizeOpts.apply=true;
F.doPrecomputation=1;
F.memoizeOpts.apply=true;

%%
Hn={F,D,C};
Fn={lklCost,mu*rglCost,cCost};
rho_n=[lklRho,rglRho,cRho];
%%
ADMM=OptiADMM([],Fn,Hn,rho_n,[]);
 maxIt=500;
 ADMM.ItUpOut=1;
% ADMM.checkConvergence=false
 ADMM.maxiter= maxIt;
ADMM.run(fft2(xinit));
%%

 