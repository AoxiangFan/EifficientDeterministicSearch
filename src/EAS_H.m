function [H, inliers] = EAS_H(X,Y,size1,size2)
D = GenerateEmbeddings(X, Y, 'A');
[~,d] = DPCP_solver(D, 2);
idx = find(d<=0.15);

th = (norm(size1) + norm(size2))*0.0016/2;
% Set option = 'LO' to activate local optimization, set option = 'None' to
% use plain RANSAC.
[H, inliers] = Post_Homography(X, Y, idx, 500, th, 'LO');