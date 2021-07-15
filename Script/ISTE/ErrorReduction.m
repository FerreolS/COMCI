%%
clear all

useGPU(1)
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
t = truth(973:972+2152,973:972+2152,:);


%% Error Reduction algorithm
Pintensity = CostComplexCircle(sz,sqrt(data)); % projection on intensity
Pphaseonly = CostComplexCircle(sz,1);
maxiter = 10000;
% init
E = zeros_(maxiter,1);
S = zeros_(maxiter,1); 
normXtrue=norm(t(:));
o = Pphaseonly.applyProx(H'*data,1);
for n=1:maxiter
    dd = H*o;
    d = Pintensity.applyProx(dd,1);
    oo = H'*d;
    o = Pphaseonly.applyProx(oo,1);
    E(n) = 20.*log10(norm(o(:)-oo(:)) +norm(d(:)-dd(:)));
    po=padarray(o,[564 564],'both');
    S(n) = 20.*log10(normXtrue/norm(t(:)-po(:)));
    disp(['iter : ', num2str(n), ' error :', num2str(E(n)), ' SNR :',num2str(S(n))])
end
