function [F, bestInliers] = Post_FundamentalMatrix(X, Y, idx, maxTrials, threshold, option)

N = size(X,1);
X = [X, ones(N,1)]';
Y = [Y, ones(N,1)]';

curTrials = 0;
bestInliers = [];
numBestInliers = 0;

while curTrials <= maxTrials
    [F, curInliers, indices] = MinimalSample2_F(X, Y, idx, threshold);
    numCurInliers = length(curInliers);
        
    if numBestInliers < numCurInliers
        bestInliers = curInliers;
        numBestInliers = numCurInliers;
        if strcmp(option,'DEGEN') || strcmp(option,'BOTH')
        % degenaracy check
            [H, degeneracy] = DegeneracyCheck(X, Y, indices, threshold);
            if degeneracy
                [curInliers, F] = DegeneracyUpdate(X, Y, threshold, bestInliers, H);
                numCurInliers = length(curInliers);
                if numBestInliers < numCurInliers
                    bestInliers = curInliers;
                    numBestInliers = numCurInliers;
                end
            end
        end
        if strcmp(option,'LO') || strcmp(option,'BOTH')
        % local optimization
            if numBestInliers >= 8
                [curInliers, F] = LocalOptimizationF(bestInliers, X, Y, threshold);
                numCurInliers = length(curInliers);
                if numBestInliers < numCurInliers
                    bestInliers = curInliers;
                    numBestInliers = numCurInliers;
                end
            end
        end
        % Update the number of trials
    end
    curTrials = curTrials + 1;
end
F = norm8Point(X(:, bestInliers), Y(:, bestInliers));
d = SampsonDistanceF(X, Y, F);
bestInliers = find(d<=threshold);
end