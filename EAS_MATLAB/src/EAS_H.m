function [H, inliers] = EAS_H(X, Y, embeddingType, size1, size2)

% remove one-to-many correspondences, which can greatly affect the
% estimation process
label = zeros(size(X,1),1);
[~, u1] = unique(Y,'rows');
X = X(u1,:);
Y = Y(u1,:);
[~, u2] = unique(X,'rows');
X = X(u2,:);
Y = Y(u2,:);

N = size(X,1);
switch embeddingType
    case 'A'
        D = GenerateEmbeddings(X, Y, embeddingType);
        [~,d] = DPCP_solver(D, 2);
        idx = find(d<=0.15);
    case 'H'
        D = GenerateEmbeddings(X, Y, embeddingType);
        [~,d] = DPCP_solver(D, 1);
        d = (d(1:2:2*N) + d(2:2:2*N))/2;
        idx = find(d<=0.15);
    case 'F'
        D = GenerateEmbeddings(X, Y, embeddingType);
        [~,d] = DPCP_solver(D, 3);
        idx = find(d<=0.15);
end

% adaptive inlier-outlier threshold in pixel, but also tunable (0 to 5 pixel in general)
th = (norm(size1) + norm(size2))*0.0016/2;

% option: local optimization (LO)
LO = true;
[H, ~] = Post_Homography(X(idx,:), Y(idx,:), 500, th, LO);
d = SampsonDistanceH([X, ones(N,1)]', [Y, ones(N,1)]', H);

label(u1(u2(d<=th))) = 1;
inliers = find(label==1);