%%
clear all
SimulateModelFieldExtension
sizeData = size(IntensityModel);


fov = sizeData(1)*dxy;
disp(['field of view (fov) width : ', num2str(fov)])
disp(['resolution at the edge of the fov : ',num2str(lambda./n0./ (fov./sqrt(fov.^2+z^2))/2)])
disp(['resolution at the center of the fov : ',num2str(lambda./n0./ (fov./2./sqrt(fov.^2./4+z^2))/2)])

disp([num2str(fov^2/(fov^2 + z^2)),'<<1 => Fresnel approximation is valid'])



%% SNR = 100
ill = 1; % Illumination
SNR = 100;
sigma = 10^(-SNR./20);
data = IntensityModel + sigma.*(random('Normal',zeros(sizeData),ones(sizeData)) + 1.i .*random('Normal',zeros(sizeData),ones(sizeData))) ;

pdata= padarray(data,[512 512],1.,'both');
sz = size(pdata);

H = LinOpPropagator(sz,lambda, n0, z,dxy,  theta, 1,1,'Fresnel');


%% Error Reduction algorithm
Pintensity = CostComplexCircle(sz,sqrt(pdata)); % projection on intensity
Pphaseonly = CostComplexCircle(sz,1);
maxiter = 100;
% init
o = Pphaseonly.applyProx(H'*pdata,1);
for n=1:maxiter
    dd = H*o;
    d = Pintensity.applyProx(dd,1);
    oo = H'*d;
    o = Pphaseonly.applyProx(oo,1);
    disp(['iter : ', num2str(n), ' error :', num2str(norm(o(:)-oo(:)) +norm(d(:)-dd(:)))])
end
% 
% %%
% 
% W = LinOpDiag(size(pdata),padarray(ones(size(data)),[512 512],0.,'both'));
% Pintensity = CostComplexCircle(sz,sqrt(pdata)); % projection on intensity
% Pphaseonly = CostComplexCircle(sz,1);
% maxiter = 100;
% % init
% o = Pphaseonly.applyProx(H'*W'*pdata,1);
% for n=1:maxiter
%     dd = W*H*o;
%     d = Pintensity.applyProx(dd,1);
%     oo = H'*W'*d;
%     o = Pphaseonly.applyProx(oo,1);
%     disp(['iter : ', num2str(n), ' error :', num2str(norm(o(:)-oo(:)) +norm(d(:)-dd(:)))])
% end
