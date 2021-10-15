function H = Hdetect(X7, Y7, F, IDXS)

ec = cross(F(:,1), F(:,2));
if (abs(ec(1)) < 1e-10) && (abs(ec(2)) < 1e-10) && (abs(ec(3)) < 1e-10)
    ec = cross(F(:,2), F(:,3));
end
ec = ec/ec(3);
Ex = skew_sym(ec);
A = inv(Ex)*F;

X3 = X7(:, int32(IDXS));
Y3 = Y7(:, int32(IDXS));

b = zeros(3,1);
for i = 1:3
    p1 = cross(Y3(i,:), A*(X3(i,:))');
    p1 = p1/p1(3);
    p2 = cross(Y3(i,:), ec);
    p2 = p2/p2(3);
    p2 = p2/norm(p2)^2;
    b(i) = p1*p2';
end

M = X3';
H = A - ec*(inv(M)*b)';

if isnan(H(1,1)) || isinf(H(1,1))
    H = eye(3);
end
