
%inim = double(imread('GrndTruth/USAF-1951_65nm_tilt_C.png'))/255.;
%inim = inim(4135:4134+20160,3857:3856+20160);
inim = double(imread('Data/Simus/GrndTruth/USAF-1951_104nm_tilt_C.png'))/255.;
inim = inim(521:520+20160,529:528+20160);
truth=1.-inim;
clear inim;
