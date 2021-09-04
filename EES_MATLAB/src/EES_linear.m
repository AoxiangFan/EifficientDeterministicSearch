function [v, inliers] = EES_linear(X, th)

N = size(X,1);
P = ones(N,1);
v = gwfit(P,X);

iter = 0;
E = 10;
t = 1;
beta = 100;
gamma = 0.9;
mbeta = 1/10;

while t > 1e-6 && iter < 100
    E0 = E;
    r = abs(X*v);
    fr = r/th;
    beta = max(gamma*beta,mbeta);
    P = exp(-fr/beta)./(exp(-fr/beta) + exp(-1/beta));
    v = gwfit(P,X);
    iter = iter + 1;
    E = P'*fr - sum(P);
    t = abs(E0 - E);
end
r = abs(X*v);
inliers = find(r<=th);