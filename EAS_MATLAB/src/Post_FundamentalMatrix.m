function [F, bestInliers] = Post_FundamentalMatrix(X, Y, maxTrials, threshold, LO, DEGEN)

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
    [F, indices, curInliers] = MinimalSample_F(X, Y, N, threshold);
    numCurInliers = length(curInliers);
    if numGoodInliers < numCurInliers
        numGoodInliers = numCurInliers;
        maxTrials = updateMaxTrials(numGoodInliers, maxTrials, N, 0.999999, 7);
        if numBestInliers < numCurInliers
            bestInliers = curInliers;
            numBestInliers = numCurInliers;
        end
        if LO
            if numCurInliers >= 16
                [F, curInliers] = LO_F(X, Y, F, threshold, 3, 50);
                numCurInliers = length(curInliers);
                if numBestInliers < numCurInliers
                    bestInliers = curInliers;
                    numBestInliers = numCurInliers;
                end
            end
        end
        if DEGEN
            [flag, H] = checkSample(X(:, indices), Y(:, indices), F, 2*threshold);
            if flag
                dH = SampsonDistanceH_all(X, Y, H);
                inliersH = find(dH <= 3*threshold);
                if length(inliersH) >= 8
                    [H, inliersH] = LO_H(X, Y, H, threshold, 3, 50);
                    if length(inliersH) >= 6
                        [F, inliersF] = H2F(X, Y, H, threshold);
                        if numBestInliers < length(inliersF)
                            bestInliers = inliersF;
                            numBestInliers = length(inliersF);
                        end
                    end
                end
            end
        end
    end
    curTrials = curTrials + 1;
end
if length(bestInliers) >= 8    
    F = norm8Point(X(:, bestInliers), Y(:, bestInliers));
    d = SampsonDistanceF(X, Y, F);
    bestInliers = find(d<=threshold);
else
    F = ones(3,3);
end
end