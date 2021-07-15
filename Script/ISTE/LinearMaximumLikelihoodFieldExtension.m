clear all
SimulateModelFieldExtension

%%
sizeData = size(AmpModel);


fov = sizeData(1)*dxy;
disp(['field of view (fov) width : ', num2str(fov)])
disp(['resolution at the edge of the fov : ',num2str(lambda./n0./ (fov./sqrt(fov.^2+z^2))/2)])
disp(['resolution at the center of the fov : ',num2str(lambda./n0./ (fov./2./sqrt(fov.^2./4+z^2))/2)])

disp([num2str(fov^2/(fov^2 + z^2)),'<<1 => Fresnel approximation is valid'])


H = LinOpPropagator(size(AmpModel),lambda, n0, z,dxy,  theta, 1,1,'Fresnel');

%% SNR = 30
ill = 1; % Illumination
SNR = 30;
sigma = 10^(-SNR./20);
data = AmpModel + sigma.*(random('Normal',zeros(sizeData),ones(sizeData)) + 1.i .*random('Normal',zeros(sizeData),ones(sizeData))) ;



o = H'*data;

figure();
plt=subplot(2,3,1);
imshow(abs(data),[])
title('|Data| SNR=100');
subplot(2,3,4)
imshow(angle(data),[])
title('angle(Data)');
subplot(2,3,3)
imshow(abs(truth(1537:1536+1024,1537:1536+1024,:)),[])
title('|truth|');
subplot(2,3,6)
imshow(angle(truth(1537:1536+1024,1537:1536+1024,:)),[])
title('angle(truth)');
subplot(2,3,2)
imshow(abs(o),[])
title('|o|');
subplot(2,3,5)
imshow(angle(o),[])
title('angle(o)');

%% SNR = 30 padding


eta = sigma./ill;
outoffov = (sqrt( 2*z*k*pi*eta.^2)  ) ./  (k *pi * eta.^2 )./dxy;
disp(['the fov extrapolation is ',num2str(outoffov), ' pixels in each direction' ])
disp(['leading to fov of width ',num2str(2*outoffov+1024), ' ( ', num2str((2*outoffov+1024)*dxy*1e3), ' mm)'])

pdata= padarray(data,[564 564],exp(1.i*mod(z*k,2*pi)),'both');
H = LinOpPropagator(size(pdata),lambda, n0, z,dxy,  theta, 1,1,'Fresnel');

o = H'*pdata;

figure();
plt=subplot(2,3,1);
imshow(abs(pdata),[])
title('|Data| padding');
subplot(2,3,4)
imshow(angle(pdata),[])
title('angle(Data)');
subplot(2,3,3)
imshow(abs(truth(973:3124,973:3124)),[])
title('|truth|');
subplot(2,3,6)
imshow(angle(truth(973:3124,973:3124)),[])
title('angle(truth)');
subplot(2,3,2)
imshow(abs(o),[])
title('|o|');
subplot(2,3,5)
imshow(angle(o),[])
title('angle(o)');

%% SNR = 30 VMLMB
pdata= padarray(data,[564 564],0.,'both');

C = LinOpCpx(size(pdata)); % utility operator to go from complex to reals pair

% Precision matrix
W = LinOpDiag(size(pdata),padarray(ones(size(data)),[564 564],0.,'both'));

%  Function definition
LS=CostL2(size(pdata),pdata,W);  % Least-Squares data term
F=LS*H*C'; %cost function

VMLMB=OptiVMLMB(F,[],[]);  
VMLMB.ItUpOut=1; 
VMLMB.maxiter=100;                  
% max number of iterations
VMLMB.OutOp=OutputOptiSNR(1,C*truth(973:3124,973:3124),1);
VMLMB.m=2;                                     % number of memorized step in hessian approximation (one step is enough for quadratic function)
VMLMB.run(C*(exp(1.i*mod(z*k,2*pi)).*ones(size(H'*pdata))));                                  % run the algorithm 
                                  % run the algorithm 
o = C'*VMLMB.xopt;
figure(); imshow(angle(o),[mod(z*k,2*pi),mod(z*k,2*pi)+1.2])
