function f = getAlgebraicF(P,D)

Ps = P.^0.5;
[~,~,v] = svd(diag(Ps)*D, 0);
f = v(:,9);

% y = -D(:,8);
% x = [D(:,1:7),D(:,9)];
% [f,~] = wfit(y,x,P);
% f = [f(1:7);1;f(8)];
% f = f/norm(f);