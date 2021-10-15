function [F, inliers] = EAS_F(X, Y, embeddingType, size1, size2)

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

% adaptive inlier-outlier threshold in pixel, but also tunable (0 to 5 pixel in general)
th = (norm(size1) + norm(size2))*0.0016/2;

% option: local optimization (LO), degenaracy updating (DEGEN)
LO = true;
DEGEN = true;
[F, ~] = Post_FundamentalMatrix(X(idx,:), Y(idx,:), 500, th, LO, DEGEN);
d = SampsonDistanceF([X, ones(N,1)]', [Y, ones(N,1)]', F);

label(u1(u2(d<=th))) = 1;
inliers = find(label==1);