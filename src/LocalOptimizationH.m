function [bestInliers,bestH] = LocalOptimizationH(curInliers,X,Y,t)
H = norm4Point(X(:, curInliers), Y(:, curInliers));
bestInliers = curInliers;
numBestInliers = length(bestInliers);
bestH = H;
irlsSteps = 10;
s = 4;

th_multiplier = 4; th_step_size = (th_multiplier*t - t)./irlsSteps;

for loirls = 0:irlsSteps
    d = SampsonDistanceH(X, Y, H);
    loind = find(d<=(th_multiplier*t - th_step_size*loirls));
    if length(loind) >= 4
        loind2 = randsample(loind, min(s*7, length(loind)));
        w = 1./(1+3*d(loind2)/t);
        H = weightedNorm4Point(X(:, loind2), Y(:, loind2), w);
        d = SampsonDistanceH(X, Y, H);
        loind = find(d<=t);
        if length(loind) > numBestInliers
            bestInliers = loind;
            numBestInliers = length(bestInliers);
            bestH = H;
        end
    end
end
end