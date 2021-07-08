%%
clear all
SimulateModel

sizeData = size(IntensityModel);


fov = sizeData(1)*dxy;
disp(['field of view (fov) width : ', num2str(fov)])
disp(['resolution at the edge of the fov : ',num2str(lambda./n0./ (fov./sqrt(fov.^2+z^2))/2)])
disp(['resolution at the center of the fov : ',num2str(lambda./n0./ (fov./2./sqrt(fov.^2./4+z^2))/2)])

disp([num2str(fov^2/(fov^2 + z^2)),'<<1 => Fresnel approximation is valid'])

H = LinOpPropagator(size(IntensityModel),lambda, n0, z,dxy,  theta, 1,1,'Fresnel');
sz = size(IntensityModel);

%% No noise no field extension
Pintensity = CostComplexCircle(sz,sqrt(IntensityModel)); % projection on intensity
Pphaseonly = CostComplexCircle(sz,1);
maxiter = 100;
% init
no = Pphaseonly.applyProx(H'*IntensityModel,1);
for n=1:maxiter
    dd = H*no;
    d = Pintensity.applyProx(dd,1);
    oo = H'*d;
    no = Pphaseonly.applyProx(oo,1);
    disp(['iter : ', num2str(n), ' error :', num2str(norm(no(:)-oo(:)) +norm(d(:)-dd(:)))])
end


%%
SimulateModelFieldExtension
%% SNR = 30
ill = 1; % Illumination
SNR = 30;
sigma = 10^(-SNR./20);
data = max(0.,IntensityModel + sigma.*(random('Normal',zeros(sizeData),ones(sizeData)))) ;
sz = size(data);
tc = truth(1537:1536+1024,1537:1536+1024,:);

%% Error Reduction algorithm
Pintensity = CostComplexCircle(sz,sqrt(data)); % projection on intensity
Pphaseonly = CostComplexCircle(sz,1);
maxiter = 10000;
% init
E = zeros_(maxiter,1);
S = zeros_(maxiter,1);
o = Pphaseonly.applyProx(H'*data,1);
for n=1:maxiter
    dd = H*o;
    d = Pintensity.applyProx(dd,1);
    oo = H'*d;
    o = Pphaseonly.applyProx(oo,1);
    E(n) = 10.*log10(norm(o(:)-oo(:)) +norm(d(:)-dd(:)));
    S(n) = 10.*log10(norm(o(:) - tc(:)));
    disp(['iter : ', num2str(n), ' error :', num2str(E(n)), ' SNR :',num2str(S(n))])
end

%% Padded Error Reduction algorithm

pdata= padarray(data,[564 564],'both');

pdatainf = padarray(sqrt(data),[564 564],0.,'both');
pdatasup = padarray(sqrt(data),[564 564],+inf,'both');

psz = size(pdata);

pH = LinOpPropagator(psz,lambda, n0, z,dxy,  theta, 1,1,'Fresnel');

%Pintensity =CostComplexCircle(psz,sqrt(pdata)); % projection on intensity
Pintensity =  CostComplexRing(psz,pdatainf,pdatasup);

Pphaseonly = CostComplexCircle(psz,1);
maxiter = 100;
% init
po = Pphaseonly.applyProx(pH'*pdata,1);

for n=1:maxiter
    dd = pH*po;
    d = Pintensity.applyProx(dd,1);
    oo = pH'*d;
    po = Pphaseonly.applyProx(oo,1);
    disp(['iter : ', num2str(n), ' error :', num2str(norm(po(:)-oo(:)) +norm(d(:)-dd(:)))])
end

%% Error Reduction algorithm with the right likelihood
Pintensity = CostIntensity(sz,sqrt(data),1./sigma.^2,'Gaussian'); % projection on intensity
Pphaseonly = CostComplexCircle(sz,1);
maxiter = 100;
% init
El = zeros(maxiter,1);
Sl = zeros(maxiter,1);
o = Pphaseonly.applyProx(H'*data,1);
for n=1:maxiter
    dd = H*o;
    d = Pintensity.applyProx(dd,1);
    oo = H'*d;
    o = Pphaseonly.applyProx(oo,1);
    El(n) = norm(o(:)-oo(:)) +norm(d(:)-dd(:));
    Sl(n) = norm(o(:) - t(:));
    disp(['iter : ', num2str(n), ' error :', num2str(norm(o(:)-oo(:)) +norm(d(:)-dd(:)))])
end
