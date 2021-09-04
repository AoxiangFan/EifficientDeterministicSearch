function [bestF, bestInlier, indices] = MinimalSample2_F(X, Y, idx, th)

indices = randsample(idx, 7);
F = norm7Point(X(:, indices), Y(:, indices));
bestF = [];
bestS = -1;
bestInlier = [];
for i = 1:size(F,3)
    d = SampsonDistanceF(X, Y, squeeze(F(:,:,i)));
    inlier = find(d<=th);
    if length(inlier) > bestS
        bestS = length(inlier);
        bestF = squeeze(F(:,:,i));
        bestInlier = inlier;
    end
end

end