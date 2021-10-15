function [bestF, bestInliers] = LO_F(X, Y, F, th, m, inlLimit)
RAN_REP = 10;

error = SampsonDistanceH(X, Y, F);
inliers = find(error<=th);
bestF = F;
bestInliers = inliers;
numBestInliers = length(bestInliers);

initialInliers = find(error<=m*th);
if length(initialInliers) <= inlLimit
    if length(inliers) >= 8
        initialF = norm8Point(X(:,initialInliers), Y(:,initialInliers));
    else
        return;
    end
else
    inliersSubset = randsample(initialInliers, inlLimit);
    initialF = norm8Point(X(:,inliersSubset), Y(:,inliersSubset));
end

error = SampsonDistanceF(X, Y, initialF);
baseInliers = find(error<=th);

if length(baseInliers) < 16
    return;
end

siz = min(int32(ceil(length(baseInliers)/2)), 14);

for i = 1:RAN_REP
    inliersSubset = randsample(baseInliers, siz);
    F = norm8Point(X(:,inliersSubset), Y(:,inliersSubset));
    [F, inliers] = iterF(X, Y, F, th, m, inlLimit);
    if length(inliers) > numBestInliers
        bestF = F;
        bestInliers = inliers;
        numBestInliers = length(inliers);
    end
end









