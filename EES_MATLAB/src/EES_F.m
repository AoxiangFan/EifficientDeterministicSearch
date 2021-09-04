function [F, inliers] = EES_F(X0, Y0, th)

N = size(X0,1);

X = [X0';ones(1,N)];
Y = [Y0';ones(1,N)];
[Xt, T1] = normalise2dpts(X);
[Yt, T2] = normalise2dpts(Y);
th = th*sqrt(T1(1,1)*T2(1,1));

D = zeros(9,N);
for nn = 1:N
    D(:,nn) = kron(Yt(:,nn),Xt(:,nn));
end
D = D';

P = ones(N,1);
f = getAlgebraicF(P,D);

iter = 0;
E = 10;
t = 1;
beta = 100;
gamma = 0.9;
mbeta = 1/10;

while t > 1e-6 && iter < 100
    E0 = E;
    r = SampsonDistanceF(Xt,Yt,reshape(f,3,3)')';
    % r = abs(C*f);
    fr = r/th;
    beta = max(gamma*beta,mbeta);
    P = exp(-fr/beta)./(exp(-fr/beta) + exp(-1.0/beta));
    f = getAlgebraicF(P,D);
    iter = iter + 1;
    E = P'*fr - sum(P);
    t = abs(E0 - E);
end
F = reshape(f,3,3)';
[u, s, v] = svd(F);
s(end) = 0;
F = u * s * v';

r = SampsonDistanceF(Xt,Yt,F)';
% r = abs(C*f);
inliers = find(r<=th);

F = T2'*F*T1;