function [Xn,m,s] = DataNorm(X)

[N, D] = size(X);
m = mean(X);
s = std(X);
Xn = X - ones(N,D)*diag(m);
Xn = Xn * (diag(1./s));

