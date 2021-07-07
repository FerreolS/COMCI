inim = double(imread('GrndTruth/USAF-1951_112nm_Smaller.png'))/255.;
%truth= exp(1i.*inim(6209:6208+23000,6441:6440+23000));
truth= exp(1i.*inim(1231:1230+23000,1866:1865+23000));

inim = [];

dxyim = 112e-9; %test image pixel size

%% Physical parameters

nAngle = 1;
IangleT =[ [  -51    -25     0     25     51         0     0     0 0];[   0     0     0     0     0     -51     -25          25 51] ];
theta = IangleT./180*pi;
lambda= 700e-9; % wavelenght
n0 = 1.52;
 z = 100*1e-6;
 k = 2*pi/lambda
dxy= 1.12e-6 ;% Camera  pixel size
 % super resolution factor 
sr=20; % dxy / dxyin
szt = size(truth);
szx = size(truth)*2;

%%
d = zeros([szx/sr,9]);
for td=1:9
    
    Iangle = theta(:,td);
    S = LinOpShape( [szx  ],[sr  szx(1)/sr  sr szx(2)/sr nAngle]);
    Frwrd = S*LinOpPropagator(szt,lambda, n0, z,dxyim,  Iangle, 2/sr,1,'BLAS2','oversample');
    
    %%
    y = Frwrd* truth;
    clear Frwrd;
    %d = squeeze( sum( sum( iy(:,223:222+256,:,223:222+256,:),3),1));
    
    d(:,:,td) = squeeze( sum( sum(abs(y(:,:,:,:,:)).^2,3),1));
    clear y;
end
save('Simu_Phase_10largeAngle_SmallZ_notilted.mat','d','truth','sr','IangleT','lambda','k','z','n0','dxy','dxyim','theta');