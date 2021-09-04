function [H, inliers] = EES_H(X0, Y0, th)

N = size(X0,1);

X = [X0';ones(1,N)];
Y = [Y0';ones(1,N)];
[Xt, T1] = normalise2dpts(X);
[Yt, T2] = normalise2dpts(Y);
th = th*sqrt(T1(1,1)*T2(1,1));

D = [];
ooo  = zeros(1,3);
for k = 1:N
  p1 = Xt(:,k);
  p2 = Yt(:,k);
  D = [ D;
    p1'*p2(3) ooo -p1'*p2(1)
    ooo p1'*p2(3) -p1'*p2(2)
   ];
end

P = ones(N,1);
H = getAlgebraicH(P,D);

iter = 0;
E = 10;
t = 1;
beta = 100;
gamma = 0.9;
mbeta = 1/10;

while t > 1e-6 && iter < 100
    E0 = E;
%     V = H*x1;
%     V = V*(diag(1./V(3,:)));
%     r = sqrt(sum((x2 - V).^2)');
    r = SampsonDistanceH(Xt,Yt,H)';
    
    fr = r/th;
    beta = max(gamma*beta,mbeta);
    P = exp(-fr/beta)./(exp(-fr/beta) + exp(-1.0/beta));
    H = getAlgebraicH(P,D);
    
    iter = iter + 1;
    E = P'*fr - sum(P);
    t = abs(E0 - E);
end

% V = H*x1;
% V = V*(diag(1./V(3,:)));
% r = sqrt(sum((x2 - V).^2)');
r = SampsonDistanceH(Xt,Yt,H)';
inliers = find(r<=th);

H = T2^(-1)*H*T1;
