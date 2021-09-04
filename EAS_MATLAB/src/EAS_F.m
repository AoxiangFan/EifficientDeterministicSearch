function [F, inliers] = EAS_F(X, Y, embeddingType, size1, size2)
N = size(X,1);
switch embeddingType
    case 'A'
        D = GenerateEmbeddings(X, Y, embeddingType);
        [~,d] = DPCP_solver(D, 2);
        idxa = find(d<=0.25);
        idxl = setdiff(1:N,idxa);
        [~,d] = DPCP_solver(D(:,idxl), 2);
        idxb = idxl(d<=0.15);
        idx = [idxa,idxb];
    case 'H'
        D = GenerateEmbeddings(X, Y, embeddingType);
        [~,d] = DPCP_solver(D, 1);
        d = (d(1:2:2*N) + d(2:2:2*N))/2;
        idxa = find(d<=0.25);
        idxl = setdiff(1:N,idxa);
        idxl2 = [];
        for i = 1:length(idxl)
            idxl2 = [idxl2,2*idxl(i)-1,2*idxl(i)];
        end
        [~,d] = DPCP_solver(D(:,idxl2), 1);
        d = (d(1:2:2*length(idxl)) + d(2:2:2*length(idxl)))/2;
        idxb = idxl(d<=0.15);
        idx = [idxa,idxb];
     case 'F'
        D = GenerateEmbeddings(X, Y, embeddingType);
        [~,d] = DPCP_solver(D, 1);
        idx = find(d<=0.15);
end

th = (norm(size1) + norm(size2))*0.0016/2;
% Set option = 'LO' to activate local optimization, option = 'DEGEN' to
% activate degeneracy check, option = 'BOTH' to activate both. Note that
% the degeneracy check part is a rough implementation and unstable, we are trying
% to improve it now. Set option  = 'None' to use plain random sample.
[F, inliers] = Post_FundamentalMatrix(X, Y, idx, 500, th, 'LO');