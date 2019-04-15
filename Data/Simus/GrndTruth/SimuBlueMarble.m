Imphase = (imread('~/matlab/COMCI/Data/Simus/GrndTruth/BlueMarble_Amp_Ph.tif','Index',2));
ImMod = (imread('~/matlab/COMCI/Data/Simus/GrndTruth/BlueMarble_Amp_Ph.tif','Index',1));

%w1 = window(@tukeywin,18000,0.002);
w1 =  1;

truth = (1.0-1./255.*double((double(ImMod).*w1) .* w1')) .* exp(1i.* pi.* double(Imphase)./180);

clear ImMod Imphase
padsz=(20160-size(truth,1))/2;
truth = padarray(truth,[padsz, padsz],'symmetric');
szt = size(truth);

% nAngle = 1;
% Iangle = zeros(2,nAngle);
% lambda= 530e-9; % wavelenght
% n0 = 1.0;
% z = 250*1e-6;
% %dxy = 5.2e-6; % Camera  pixel size
% dxy= 1.04e-6 ;
% 
% 
% % super resolution factor
% sr=10; % dxy / dxyin
% 
% S = LinOpShape( [szx  ],[sr  szx(1)/sr  sr szx(2)/sr nAngle]);
% %Scaling(1./sr)
% Frwrd =LinOpDiag([sr  szx(1)/sr  sr szx(2)/sr nAngle],1./sr)* S*LinOpPropagator(szx,lambda, n0, z,dxy/sr,  Iangle, 1,1,'AS');
% 
% y = Frwrd* data;
% iy =(abs(y)).^2;
% d1 = squeeze( sum( sum( iy(:,:,:,:,:),3),1));