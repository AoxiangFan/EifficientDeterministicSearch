function [H, inliers] = EAS_H(X, Y, embeddingType, size1, size2)
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

th = (norm(size1) + norm(size2))*0.0016/2;
% Set option = 'LO' to activate local optimization, set option = 'None' to
% use plain random sample.
[H, inliers] = Post_Homography(X, Y, idx, 500, th, 'LO');