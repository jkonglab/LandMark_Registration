clear;
addpath(genpath('.'))

ref_rgb=imread('ref.tif');
target_rgb=imread('target.tif');
ref=rgb2gray(ref_rgb);
target=rgb2gray(target_rgb);

%registration
[D,registered] = imregdemons(target,ref);
registered_rgb = imwarp(target_rgb,D);

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