function f = gwfit(P,D)

N = size(D,1);
Ps = P.^0.5;
dP = sparse(1:N,1:N,Ps,N,N,N);
[~,~,v] = svd(dP*D, 0);
f = v(:,end);

% y = -D(:,8);
% x = [D(:,1:7),D(:,9)];
% [f,~] = wfit(y,x,P);
% f = [f(1:7);1;f(8)];
% f = f/norm(f);