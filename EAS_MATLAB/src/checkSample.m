function [flag, H] = checkSample(X7, Y7, F, th)
IDXS = [1,2,3;4,5,6;1,2,7;4,5,7;3,6,7];
for i = 1:5
    H = Hdetect(X7, Y7, F, IDXS(i,:));
    d = SampsonDistanceH(X7, Y7, H);
    [~, idx] = sort(d);
    H = norm4Point(X7(:,idx(1:5)), Y7(:,idx(1:5)));
    d2 = SampsonDistanceH(X7, Y7, H);
    inliers = find(d2<=th);
    if length(inliers) > 4
        flag = true;
        return;
    end
end
flag = false;
        
