inim = 1-double(imread('GrndTruth/USAF-1951_6700nm_tilt_C.png'))/255.;
truth= padarray( padarray(exp(1i.*inim), [612+512,620+512],'symmetric','pre'), [612+512,619+512],'symmetric','post');
%truth= padarray( padarray(exp(1i.*inim), [612+512,620+512],1.,'pre'), [612+512,619+512],1.,'post');
truth= padarray(truth,[512 512],1.,'both'); 
%truth= exp(1i.*inim);
inim = [];


nAngle = 1;
Iangle =[ 0 0]';
theta = Iangle./180*pi;
lambda= 532e-9; % wavelenght
n0 = 1.0;
k = 2.*pi/lambda*n0; % wavenumber

z = 265e-3;
 
dxy = 6.7e-6; % Camera  pixel size
sr=2; % dxy / dxyin
dxyim = dxy/sr; %65e-9; %test image pixel size


 % super resolution factor 
szt = size(truth);
szx = size(truth)*2;


%%
% useGPU(1);
% gpuCpuConverter(truth);
%%
S = 1./(sr)*LinOpShape( [szx  ],[sr  szx(1)/sr  sr szx(2)/sr nAngle]);
Frwrd = S*LinOpPropagator(szt,lambda, n0, z,dxyim*2,  theta, 1,1,'BLAS','oversample');

y = Frwrd* truth;
clear Frwrd;
py = y;
pd = squeeze( sum( sum( py(:,:,:,:,:),3),1));
pd1 = pd(1537:1536+1024,1537:1536+1024,:);


y =(abs(y)).^2;

d = squeeze( sum( sum( y(:,:,:,:,:),3),1));
d1 = d(1537:1536+1024,1537:1536+1024,:);