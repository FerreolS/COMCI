inim = double(imread('GrndTruth/USAF-1951_140nm_tilt_C.png'))/255.;
truth= exp(1i.*inim(5731:5730+23000,4971:4970+23000));
inim = [];
dxyim = 140e-9; %test image pixel size


nAngle = 1;
IangleT =[ [  -50    -25     0     25     50         0     0     0 0];[   0     0     0     0     0     -50     -25          25 50] ];
IangleT = IangleT./180*pi;
lambda= 700e-9; % wavelenght
n0 = 1.0;
 z = 100*1e-6;
 
dxy = 5.2e-6; % Camera  pixel size
dxy= 1.12e-6 ;
 % super resolution factor 
sr=16; % dxy / dxyin
szt = size(truth);
szx = size(truth)*2;

%%
d = zeros([szx/sr,9]);
for td=1:9
    
    Iangle = IangleT(:,td);
    S = LinOpShape( [szx  ],[sr  szx(1)/sr  sr szx(2)/sr nAngle]);
    Frwrd = S*LinOpPropagator(szt,lambda, n0, z,dxyim,  Iangle, 1,1,'BLAS2','oversample');
    
    %%
    y = Frwrd* truth;
    clear Frwrd;
    y =(abs(y)).^2;
    %d = squeeze( sum( sum( iy(:,223:222+256,:,223:222+256,:),3),1));
    
    d(:,:,td) = squeeze( sum( sum( y(:,:,:,:,:),3),1));
    clear y;
end
%d1 = d(925:924+1024,925:924+1024,:);
%save('Angle_AS_Simu.mat','d','d1','sr','IangleT','lambda','k','z','n0','dxy','dxyim');
