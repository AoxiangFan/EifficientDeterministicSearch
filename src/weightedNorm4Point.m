function H = weightedNorm4Point(X,Y,w)
[~,c] = size(X);

if (size(X,1) == 2)
  X = [X ; ones(1,size(X,2))];
  Y = [Y ; ones(1,size(Y,2))];
end

[X, T1] = normalise2dpts(X);
[Y, T2] = normalise2dpts(Y);
X(isnan(X)) = 1;
Y(isnan(Y)) = 1;

D = [];
ooo  = zeros(1,3);
for k=1:c
  p1 = X(:,k);
  p2 = Y(:,k);
  D = [ D;
    p1'*p2(3) ooo -p1'*p2(1)
    ooo p1'*p2(3) -p1'*p2(2)
   ];
end
w = kron(diag(sqrt(w)),[1 0;0 1]);
[~,~,v] = svd(w*D, 0);
h = v(:,9);
H = reshape(h,3,3)';

H = T2^(-1)*H*T1;
end