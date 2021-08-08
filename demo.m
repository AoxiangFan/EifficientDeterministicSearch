clear,clc
close all

%% example for fundamental matrix estimation

% A normal example
load ./data/wash.mat

% A degeneracy example
% load ./data/box.mat 

size1 = [size(I1,1), size(I1,2)];
size2 = [size(I2,1), size(I2,2)];

[F, inliers] = EAS_F(X0, Y0, size1, size2);
mean_error = mean(SampsonDistanceF(vpts(1:3,:),vpts(4:6,:),F))
figure
showMatchedFeatures(I1, I2, X0(inliers,:), Y0(inliers,:), 'montage');

%% example for homography estimation

load ./data/graf.mat
size1 = [size(I1,1), size(I1,2)];
size2 = [size(I2,1), size(I2,2)];

[H, inliers] = EAS_H(X0, Y0, size1, size2);
mean_error = mean(SampsonDistanceH(vpts(1:3,:),vpts(4:6,:),H))
figure
showMatchedFeatures(I1, I2, X0(inliers,:), Y0(inliers,:), 'montage');
