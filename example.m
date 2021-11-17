clear;
addpath(genpath('.'))
ref_rgb=imread('ref.tif');
target_rgb=imread('target.tif');
eta=0.25;lr=200;lt=500;tau=0.1;n=4; m=2;nuclei=1;

% registration
[registered_rgb,delta_x,delta_y]=register(ref_rgb,target_rgb,eta,lr,lt,tau,n,m,nuclei);

ref=rgb2gray(ref_rgb);
target=rgb2gray(target_rgb);
registered=rgb2gray(registered_rgb);

% metrics before registration
ssim1=ssim(ref,target);
mse1=immse(ref,target);
mi1=MutualInformation(ref,target);
cor1=corr2(ref,target);

% metrics after registration
ssim2=ssim(ref,registered);
mse2=immse(ref,registered);
mi2=MutualInformation(ref,registered);
cor2=corr2(ref,registered);