function [bestF, bestInliers] = innerFH(Xh, Yh, Xo, Yo, X, Y, th, num, sam_sizH, sam_sizO)
lenH = size(Xh,2);
lenO = size(Xo,2);
max_i = 0;
max_s = 0;

bestF = [];
bestInliers = [];

for i = 1:num
    hsam = randperm(lenH, sam_sizH);
    osam = randperm(lenO, sam_sizO);
    Xh_sam = Xh(:,hsam);
    Yh_sam = Yh(:,hsam);
    Xo_sam = Xo(:,osam);
    Yo_sam = Yo(:,osam);
    
    X_sam = [Xh_sam, Xo_sam];
    Y_sam = [Yh_sam, Yo_sam];
    
    aF = norm8Point(X_sam, Y_sam);
    d = SampsonDistanceF(X, Y, aF);
    inliers = find(d<=th);
    nI = length(inliers);
    
    if max_i < nI
        bestInliers = inliers;
        bestF = aF;
        max_i = nI;
    end
    
    if nI > max_s
        max_s = nI;
        [aF, inliers] = iterF(X, Y, aF, th, 2, 50);
        nI = length(inliers);
    end
    
    if max_i < nI
        bestInliers = inliers;
        bestF = aF;
        max_i = nI;
    end
end
        
        
        
        
        
        
        