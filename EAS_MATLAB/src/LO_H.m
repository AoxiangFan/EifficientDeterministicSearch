function [bestH, bestInliers] = LO_H(X, Y, H, th, m, inlLimit)
RAN_REP = 10;

error = SampsonDistanceH_all(X, Y, H);
inliers = find(error<=th);
bestH = H;
bestInliers = inliers;
numBestInliers = length(bestInliers);

initialInliers = find(error<=m*th);
if length(initialInliers) <= inlLimit
    if length(inliers) >= 8
        initialH = norm4Point(X(:,initialInliers), Y(:,initialInliers));
    else
        return;
    end
else
    inliersSubset = randsample(initialInliers, inlLimit);
    initialH = norm4Point(X(:,inliersSubset), Y(:,inliersSubset));
end


error = SampsonDistanceH_all(X, Y, initialH);
baseInliers = find(error<=th);

if length(baseInliers) < 16
    return;
end

siz = min(int32(ceil(length(baseInliers)/2)), 14);

for i = 1:RAN_REP
    inliersSubset = randsample(baseInliers, siz);
    H = norm4Point(X(:,inliersSubset), Y(:,inliersSubset));
    [H, inliers] = iterH(X, Y, H, th, m, inlLimit);
    if length(inliers) > numBestInliers
        bestH = H;
        bestInliers = inliers;
        numBestInliers = length(inliers);
    end
end









