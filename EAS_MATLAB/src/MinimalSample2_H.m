function [H, Inlier, indices] = MinimalSample2_H(X, Y, idx, th)

indices = randsample(idx, 4);
H = norm4Point(X(:, indices), Y(:, indices));
d = SampsonDistanceH(X, Y, H);
Inlier = find(d<=th);

end