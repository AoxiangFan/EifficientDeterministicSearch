clear,clc
close all

%% example for synthetic data

disp('synthetic example: ');

% experiment settings
outlier_rate = 0.8;
D = 8; d = 7; % dimension of data
N = 500; % number of inliers
ratio = 1 ./ (1 ./ outlier_rate - 1);
M = floor(N * ratio);
sigma = 0.15; % noise level for inliers

% generate synthetic data
S = orth(randn(D,d));
X = normc( S*randn(d,N) + sigma * randn(D,N) );
O = normc(randn(D,M));
Xtilde = [X O];
vg = gwfit(ones(d,1),S'); % ground truth vector

% run EES
th = 0.8*sigma;
[ve1, ~] = EES_linear(Xtilde', th);

% run DPCP method from NeurIPS2018 paper "Dual principal component pursuit: Improved analysis and efficient algorithms"
[ve2, ~] = DPCP_solver(Xtilde,1);

% angular error
disp(['EES angular error: ',num2str(acos(norm(vg'*ve1))*180/pi),' degrees']);
disp(['DPCP angular error: ',num2str(acos(norm(vg'*ve2))*180/pi),' degrees']);

%% example for road plane estimation
% The original road plane data and visualization code are from the release of
% NeurIPS2018 paper "Dual principal component pursuit: Improved analysis and efficient algorithms".

disp('road plane estimation example: ');

% load and prepare data
load('./data/005_0153_data.mat');

sigma = 0.08; % noise level
inliers = inliers + [sigma*randn(3,size(inliers,2));zeros(1,size(inliers,2))];
outliers = outliers + [sigma*randn(3,size(outliers,2));zeros(1,size(outliers,2))];

X_tilde = [inliers outliers];

% run EES
[ve1, ~] = EES_linear(X_tilde', gt_threshold);
[inliers1, outliers1] = cut(X_tilde, ve1, gt_threshold);
inliers1 = X_tilde(:,inliers1);
outliers1 = X_tilde(:,outliers1);

% run DPCP
[ve2, ~] = DPCP_solver(X_tilde,1);

[inliers2, outliers2] = cut(X_tilde, ve2, gt_threshold);
inliers2 = X_tilde(:,inliers2);
outliers2 = X_tilde(:,outliers2);


% visualization
i = 153;
img = imread('./data/0000000153.png');

figure;
imshow(Draw_image(img,inliers,outliers));
text(20,40,'Ground Truth','FontSize',50,'Color','white');

figure;
imshow(Draw_image(img,inliers2,outliers2));
text(20,40,'DPCP','FontSize',50,'Color','white');

figure;
imshow(Draw_image(img,inliers1,outliers1));
text(20,40,'EES','FontSize',50,'Color','white');







function [Roadimage] = Draw_image(img,inliers,outliers)
    %% Camera matrixs
    velo2rect = [ 2.34773698e-04,1.04494074e-02,9.99945389e-01,0.00000000e+00;
                 -9.99944155e-01,1.05653536e-02,1.24365378e-04,0.00000000e+00;
                 -1.05634778e-02,-9.99889574e-01,1.04513030e-02,0.00000000e+00;
                 -2.79681694e-03,-7.51087914e-02,-2.72132796e-01,1.00000000e+00];

    rect2disp = [721.5377,0.,0.,0.;
                   0.,721.5377,0.,0.;
                 609.5593 172.854,0.,1.;
                   0.,0.,387.5744,0.];
    %% inliers
    inliers_img = inliers'*velo2rect*rect2disp;
    inliers_img = inliers_img';
    inliers_img = inliers_img./inliers_img(4,:);

    inliers_visiable_idx = find(inliers_img(1,:) <= 1242 &...
                                    inliers_img(1,:) >= 0 &...
                                    inliers_img(2,:) <= 375 &...
                                    inliers_img(2,:) >= 0&...
                                    inliers_img(3,:) <= 255&...
                                    inliers_img(3,:) >= 0);
    inliers_img = inliers_img(:,inliers_visiable_idx);

    inliers2d_img = round(inliers_img(1:2,:))+1;
    for i = 1:size(inliers2d_img,2)
        y = inliers2d_img(1,i);
        x = inliers2d_img(2,i);
        img(x,y,1) = 0;
        img(x,y,2) = 0;
        img(x,y,3) = 255;
    end
    %% outliers
    outliers_img = outliers'*velo2rect*rect2disp;
    outliers_img = outliers_img';
    outliers_img = outliers_img./outliers_img(4,:);

    outliers_visiable_idx = find(outliers_img(1,:) <= 1242 &...
                                    outliers_img(1,:) >= 0 &...
                                    outliers_img(2,:) <= 375 &...
                                    outliers_img(2,:) >= 0&...
                                    outliers_img(3,:) <= 255&...
                                    outliers_img(3,:) >= 0);
    outliers_img = outliers_img(:,outliers_visiable_idx);

    outliers2d_img = round(outliers_img(1:2,:))+1;
    for i = 1:size(outliers2d_img,2)
        y = outliers2d_img(1,i);
        x = outliers2d_img(2,i);
        img(x,y,1) = 255;
        img(x,y,2) = 0;
        img(x,y,3) = 0;
    end
    Roadimage = img;
end

function [idInliers, idOutliers] = cut(X_tilde, B, th)
    dists = abs(B'*normc(X_tilde));
    idInliers = find(dists <= th);
    idOutliers = find(dists > th);
end


