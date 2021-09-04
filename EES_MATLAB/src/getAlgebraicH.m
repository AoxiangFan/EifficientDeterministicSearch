function H = getAlgebraicH(P,D)

Ps = repmat(P,1,2)';
w = Ps(:);
y = D(:,end);
x = D(:,1:end-1);
[h,~] = wfit(y,x,w);
h = [-h;1];
% [~,~,v] = svd(W*D, 0);
% h = v(:,9);

H = reshape(h,3,3)';
% a = W*D*h;
% 
% N = length(P);
% r2 = zeros(N,1);
% for i = 1:N
%     r2(i) = a(2*i-1)^2 + a(2*i)^2;
% end