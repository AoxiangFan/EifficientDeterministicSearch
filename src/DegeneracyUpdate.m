function [bestInliers, bestF] = DegeneracyUpdate(X, Y, th, bestInliers, H)

th_t  = max(th, 5);
d = SampsonDistanceH(X,Y,H);
ind_a = find(d>2.0*th_t);
ind_b = find(d<=th_t);
bestF = rand(3,3);
if length(ind_a) >= 2 && length(ind_b) >= 6
    bestS = 0;
    for i = 1:100
        sp_a = randsample(ind_a,2);        
        sp_b = randsample(ind_b,6);
        sample = [sp_a,sp_b];
        F = norm8Point(X(:,sample),Y(:,sample));
        d = SampsonDistanceF(X,Y,F);
        curS = length(find(d<=th_t));
        if curS > bestS
            bestS = curS;
            bestInliers = find(d<=th_t);
            bestF = F;
        end
    end
end

