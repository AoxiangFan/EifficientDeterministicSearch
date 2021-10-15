function [H, indices] = MinimalSample_H(X, Y, N)

indices = randperm(N, 4);
H = norm4Point(X(:, indices), Y(:, indices));
end