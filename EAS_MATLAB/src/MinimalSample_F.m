function [bestF, indices, bestInliers] = MinimalSample_F(X, Y, N, th)

indices = randperm(N, 7);
F = norm7Point(X(:, indices), Y(:, indices));
bestF = [];
bestS = -1;
bestInliers = [];
for i = 1:size(F,3)
    d = SampsonDistanceF(X, Y, squeeze(F(:,:,i)));
    inliers = find(d<=th);
    if length(inliers) > bestS
        bestS = length(inliers);
        bestF = squeeze(F(:,:,i));
        bestInliers = inliers;
    end
end

end