function [F, inliers] = EAS_F(X,Y,size1,size2)
D = GenerateEmbeddings(X, Y, 'A');
N = size(D,2);
[~,d] = DPCP_solver(D, 2);
idxa = find(d<=0.25);
idxl = setdiff(1:N,idxa);
[~,d] = DPCP_solver(D(:,idxl), 2);
idxb = idxl(d<=0.15);
idx = [idxa,idxb];

th = (norm(size1) + norm(size2))*0.0016/2;
% Set option = 'LO' to activate local optimization, option = 'DEGEN' to
% activate degeneracy check, option = 'BOTH' to activate both. Note that
% the degeneracy check part is a rough implementation and unstable, we are trying
% to improve it now. Set option  = 'None' to use plain RANSAC.
[F, inliers] = Post_FundamentalMatrix(X, Y, idx, 500, th, 'LO');