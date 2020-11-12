load('TestSize3_x100.mat')

figure(12);hold off
figure(12);plot(wdSz*dxy*1e6,squeeze(P1(:,1,11)/nb+20))
hold on
figure(12);plot(wdSz*dxy*1e6,squeeze(P1(:,2,11)/nb+20))
figure(12);plot(wdSz*dxy*1e6,squeeze(P1(:,3,11)/nb+20))
figure(12);plot(wdSz*dxy*1e6,squeeze(P1(:,5,11)/nb+20))
theta = IncidenceA *pi/180
 SNR11 = SNR(11)

 Bnd = (sqrt( 2*z*k*pi*SNR11.^2*cos(theta)+ sin(theta).^2) + sin(abs(theta)) ) ./  (k *pi * SNR11.^2 * cos(theta).^2).*1e6
 axis([0 250 0 5])
line([Bnd(1) Bnd(1)],[0 1.5],'Color',[1 0 0])
line([Bnd(2) Bnd(2)],[0 1],'Color',[1 0 0])
line([Bnd(3) Bnd(3)],[0 1],'Color',[1 0 0])
line([Bnd(5) Bnd(5)],[0 1],'Color',[1 0 0])

%%
figure(13);hold off
figure(13);plot(wdSz*dxy*1e6,squeeze(P1(:,2,9))/nb-20*log10(SNR(9)));
hold on
figure(13);plot(wdSz*dxy*1e6,squeeze(P1(:,2,10))/nb-20*log10(SNR(10)));
figure(13);plot(wdSz*dxy*1e6,squeeze(P1(:,2,11))/nb-20*log10(SNR(11)));
figure(13);plot(wdSz*dxy*1e6,squeeze(P1(:,2,12))/nb-20*log10(SNR(12)));
figure(13);plot(wdSz*dxy*1e6,squeeze(P1(:,2,13))/nb-20*log10(SNR(13)));

theta45 = IncidenceA(2) *pi/180
SNR = SNR(9:13);
 Bnd = (sqrt( 2*z*k*pi.*SNR.^2.*cos(theta45)+ sin(theta45).^2) + sin(abs(theta45)) ) ./  (k *pi .* SNR.^2 .* cos(theta45).^2).*1e6
axis([0 250 0 5])
line([Bnd(1) Bnd(1)],[0 1.5],'Color',[1 0 0])
line([Bnd(2) Bnd(2)],[0 1],'Color',[1 0 0])
line([Bnd(3) Bnd(3)],[0 1],'Color',[1 0 0])
line([Bnd(4) Bnd(4)],[0 1],'Color',[1 0 0])
line([Bnd(5) Bnd(5)],[0 1],'Color',[1 0 0])
%%
load('TestSpectralCoherence_Gauss30.mat')

Bndc = CoherenceLength ./ cos(-30/180*pi)^2.*(sqrt(1+2.*cos(-30/180*pi).*z ./CoherenceLength) + sin(-30/180*pi))*1e6


hold off
figure(25);plot(wdSz*1e6*dxy,squeeze(PSpectCoherence(:,7,1)/nb-PSpectCoherence(end,7,1)/nb))
axis([0 250 -0.6 9])
hold on
%figure(25);plot(wdSz*1e6*dxy,squeeze(PSpectCoherence(:,6,1)/nb-PSpectCoherence(end,6,1)/nb))
figure(25);plot(wdSz*1e6*dxy,squeeze(PSpectCoherence(:,5,1)/nb-PSpectCoherence(end,5,1)/nb))
figure(25);plot(wdSz*1e6*dxy,squeeze(PSpectCoherence(:,1,1)/nb-PSpectCoherence(end,1,1)/nb))
line([Bndc(end) Bndc(end)],[-0.5 0.1],'Color',[1 0 0])
%line([Bndc(6) Bndc(6)],[-0.5 0.5],'Color',[1 0 0])
line([Bndc(5) Bndc(5)],[-0.5 0.5],'Color',[1 0 0])
line([Bndc(1) Bndc(1)],[-0.5 0.5],'Color',[1 0 0])


%%
load('TestSpatialCoherence_Disk_Angle30.mat')
figure(26);plot(wdSz*1e6*dxy,squeeze(PSpatCoherence(:,1,1)/nb-PSpatCoherence(end,1,1)/nb))
axis([0 250 -0.8 8])
hold on
figure(26);plot(wdSz*1e6*dxy,squeeze(PSpatCoherence(:,2,1)/nb-PSpatCoherence(end,2,1)/nb))
figure(26);plot(wdSz*1e6*dxy,squeeze(PSpatCoherence(:,3,1)/nb-PSpatCoherence(end,3,1)/nb))
figure(26);plot(wdSz*1e6*dxy,squeeze(PSpatCoherence(:,4,1)/nb-PSpatCoherence(end,4,1)/nb))

Bndc = lambda ./ tan( HalfApexAngle)./cos(IncidenceA /180*pi)*1e6/pi;
line([Bndc(end) Bndc(end)],[-0.8 0.1],'Color',[1 0 0])
line([Bndc(3) Bndc(3)],[-0.8 0.5],'Color',[1 0 0])
line([Bndc(2) Bndc(2)],[-0.8 0.5],'Color',[1 0 0])
line([Bndc(1) Bndc(1)],[-0.8 0.5],'Color',[1 0 0])



hold off;


%%
load('TestSRCohe4x100.mat')
figure(28);hold on;
axis([0 8 -43 -32])
t = 1;plot(dxy.*sr*1e6,Psr1(:,t)/nb)
t = 2;plot(dxy.*sr*1e6,Psr1(:,t)/nb)
Bndc = CoherenceLength ./ cos(-0/180*pi)^2.*(sqrt(1+2.*cos(-0/180*pi).*z ./CoherenceLength) + sin(-0/180*pi))
pixsz = lambda./ (Bndc./sqrt(Bndc.^2+z^2))*1e6/2
line([pixsz(2) pixsz(2)],[-42 -40],'Color',[1 0 0])
line([pixsz(1) pixsz(1)],[-42 -40],'Color',[1 0 0])


hold off;
