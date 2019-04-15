
inim = double(imread('GrndTruth/USAF-1951_65nm_tilt_C.png'))/255.;
inim = inim(4135:4134+20160,3857:3856+20160);
truth=1.-inim;
clear inim;
