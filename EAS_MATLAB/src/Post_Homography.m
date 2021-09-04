function [H, bestInliers] = Post_Homography(X, Y, idx, maxTrials, threshold, option)

N = size(X,1);
X = [X, ones(N,1)]';
Y = [Y, ones(N,1)]';

curTrials = 0;
bestInliers = [];
numBestInliers = 0;

while curTrials <= maxTrials
    [H, curInliers, ~] = MinimalSample2_H(X, Y, idx, threshold);
    numCurInliers = length(curInliers);
        
    if numBestInliers < numCurInliers
        bestInliers = curInliers;
        numBestInliers = numCurInliers;
        if strcmp(option, 'LO')
        % local optimization
            if numBestInliers >= 4
                [curInliers, H] = LocalOptimizationH(bestInliers, X, Y, threshold);
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
H = norm4Point(X(:, bestInliers), Y(:, bestInliers));
d = SampsonDistanceH(X, Y, H);
bestInliers = find(d<=threshold);
end