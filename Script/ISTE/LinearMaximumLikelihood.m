clear all
SimulateModel

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
title('|Data| SNR=40');
subplot(2,3,4)
imshow(angle(data),[])
title('angle(Data)');
subplot(2,3,3)
imshow(abs(truth),[])
title('|truth|');
subplot(2,3,6)
imshow(angle(truth),[])
title('angle(truth)');
subplot(2,3,2)
imshow(abs(o),[])
title('|o|');
subplot(2,3,5)
imshow(angle(o),[])
title('angle(o)');
%% SNR = 20 
SNR = 20;
sigma = 10^(-SNR./20);
data = AmpModel + sigma.*(random('Normal',zeros(sizeData),ones(sizeData)) + 1.i .*random('Normal',zeros(sizeData),ones(sizeData))) ;


o = H'*data;

figure();
plt=subplot(2,3,1);
imshow(abs(data),[])
title('|Data| SNR=20');
subplot(2,3,4)
imshow(angle(data),[])
title('angle(Data)');
subplot(2,3,3)
imshow(abs(truth),[])
title('|truth|');
subplot(2,3,6)
imshow(angle(truth),[])
title('angle(truth)');
subplot(2,3,2)
imshow(abs(o),[])
title('|o|');
subplot(2,3,5)
imshow(angle(o),[])
title('angle(o)');

%% SNR = 10
SNR=10;
sigma = 10^(-SNR./20);
data = AmpModel + sigma.*(random('Normal',zeros(sizeData),ones(sizeData)) + 1.i .*random('Normal',zeros(sizeData),ones(sizeData))) ;


o = H'*data;

figure();
plt=subplot(2,3,1);
imshow(abs(data),[])
title('|Data| SNR=1');
subplot(2,3,4)
imshow(angle(data),[])
title('angle(Data)');
subplot(2,3,3)
imshow(abs(truth),[])
title('|truth|');
subplot(2,3,6)
imshow(angle(truth),[])
title('angle(truth)');
subplot(2,3,2)
imshow(abs(o),[])
title('|o|');
subplot(2,3,5)
imshow(angle(o),[])
title('angle(o)');
