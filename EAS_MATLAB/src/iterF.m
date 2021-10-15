function [bestF, bestInliers] = iterF(X, Y, F, th, m, inlLimit)
ILSQ_ITERS = 4;
ths = m*th;
dth = (ths - th)/ILSQ_ITERS;

error = SampsonDistanceF(X, Y, F);
inliers = find(error<=th);

numBestInliers = length(inliers);
bestInliers = inliers;
bestF = F;

if length(inliers) <= inlLimit
    if length(inliers) >= 8
        F = norm8Point(X(:,inliers), Y(:,inliers));
    else
        return;
    end
else
    inliersSubset = randsample(inliers, inlLimit);
    F = norm8Point(X(:,inliersSubset), Y(:,inliersSubset));
end

for i = 1:ILSQ_ITERS
    error = SampsonDistanceF(X, Y, F);
    inliersA = find(error<=th);
    inliersB = find(error<=ths);
    if length(inliersA) > numBestInliers
        numBestInliers = length(inliersA);
        bestF = F;
        bestInliers = inliersA;
    end
    if length(inliersB) <= inlLimit
        if length(inliersB) >= 8
            F = norm8Point(X(:, inliersB), Y(:, inliersB));
        else
            return;
        end
    else
        inliersSubset = randsample(inliersB, inlLimit);
        F = norm8Point(X(:,inliersSubset), Y(:,inliersSubset));    
    end
    ths = ths - dth;
end
    
    
    
    
    
    
    
    
    
        

