function H = norm4Point(xs1,xs2)
[~,c] = size(xs1);

if (size(xs1) ~= size(xs2))
 error ('Input point sets are different sizes!')
end

if (size(xs1,1) == 2)
  xs1 = [xs1 ; ones(1,size(xs1,2))];
  xs2 = [xs2 ; ones(1,size(xs2,2))];
end

[xs1, T1] = normalise2dpts(xs1);
[xs2, T2] = normalise2dpts(xs2);
xs1(isnan(xs1)) = 1;
xs2(isnan(xs2)) = 1;

D = [];
ooo  = zeros(1,3);
for k=1:c
  p1 = xs1(:,k);
  p2 = xs2(:,k);
  D = [ D;
    p1'*p2(3) ooo -p1'*p2(1)
    ooo p1'*p2(3) -p1'*p2(2)
   ];
end

[~,~,v] = svd(D, 0);
h = v(:,9);
H = reshape(h,3,3)';

H = T2^(-1)*H*T1;

end