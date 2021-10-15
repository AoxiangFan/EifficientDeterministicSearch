function [bestH, bestInliers] = iterH(X, Y, H, th, m, inlLimit)
ILSQ_ITERS = 4;
ths = m*th;
dth = (ths - th)/ILSQ_ITERS;

error = SampsonDistanceH_all(X, Y, H);
inliers = find(error<=th);

numBestInliers = length(inliers);
bestInliers = inliers;
bestH = H;

if length(inliers) <= inlLimit
    if length(inliers) >= 4
        H = norm4Point(X(:,inliers), Y(:,inliers));
    else
        return;
    end
else
    inliersSubset = randsample(inliers, inlLimit);
    H = norm4Point(X(:,inliersSubset), Y(:,inliersSubset));
end

for i = 1:ILSQ_ITERS
    error = SampsonDistanceH_all(X, Y, H);
    inliersA = find(error<=th);
    inliersB = find(error<=ths);
    if length(inliersA) > numBestInliers
        numBestInliers = length(inliersA);
        bestH = H;
        bestInliers = inliersA;
    end
    if length(inliersB) <= inlLimit
        if length(inliersB) >= 4
            H = norm4Point(X(:, inliersB), Y(:, inliersB));
        else
            return;
        end
    else
        inliersSubset = randsample(inliersB, inlLimit);
        H = norm4Point(X(:,inliersSubset), Y(:,inliersSubset));    
    end
    ths = ths - dth;
end
    
    
    
    
    
    
    
    
    
        

