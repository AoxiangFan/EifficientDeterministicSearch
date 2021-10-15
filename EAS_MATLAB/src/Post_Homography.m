function [H, bestInliers] = Post_Homography(X, Y, maxTrials, threshold, LO)

N = size(X,1);
X = [X, ones(N,1)]';
Y = [Y, ones(N,1)]';

curTrials = 0;
bestInliers = [];
numBestInliers = 0;

numGoodInliers = 0;

global Dx;
global Dy;
[Dx, Dy] = preSampsonDistanceH_all(X, Y);

while curTrials <= maxTrials
    [H, ~] = MinimalSample_H(X, Y, N);
    d = SampsonDistanceH_all(X, Y, H);
    curInliers = find(d<=threshold);
    numCurInliers = length(curInliers);
    if numGoodInliers < numCurInliers
        numGoodInliers = numCurInliers;
        maxTrials = updateMaxTrials(numGoodInliers, maxTrials, N, 0.999999, 4);
        if numBestInliers < numCurInliers
            bestInliers = curInliers;
            numBestInliers = numCurInliers;
        end
        if LO
            if numCurInliers >= 8
                [H, curInliers] = LO_H(X, Y, H, threshold, 3, 50);
                numCurInliers = length(curInliers);
                if numBestInliers < numCurInliers
                    bestInliers = curInliers;
                    numBestInliers = numCurInliers;
                end
            end
        end
    end
    curTrials = curTrials + 1;
end
if length(bestInliers) >= 4    
    H = norm4Point(X(:, bestInliers), Y(:, bestInliers));
    d = SampsonDistanceH_all(X, Y, H);
    bestInliers = find(d<=threshold);
else
    H = ones(3,3);
end
end