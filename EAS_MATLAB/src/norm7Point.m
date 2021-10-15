function Fm = norm7Point(X, Y)

num = size(X, 2);
[X, T1] = normalise2dpts(X);
[Y, T2] = normalise2dpts(Y);
if any(any(isnan(X))) || any(any(isnan(Y)))
    Fm = rand(3,3);
    return
end

m = zeros(num, 9);
for idx = 1: num
  m(idx,:) = [...
    X(1,idx)*Y(1,idx), X(2,idx)*Y(1,idx), Y(1,idx), ...
    X(1,idx)*Y(2,idx), X(2,idx)*Y(2,idx), Y(2,idx), ...
                 X(1,idx),              X(2,idx), 1];
end
% m(:,1:3) = (repmat(pts2h(1,:), 3, 1).*pts1h)';
% m(:,4:6) = (repmat(pts2h(2,:), 3, 1).*pts1h)';
% m(:,7:9) = pts1h';


[~, ~, vm] = svd(m, 0);

F = reshape(vm(:, end), 3, 3)';
f11=F(1,1);f12=F(1,2);f13=F(1,3);f21=F(2,1);f22=F(2,2);f23=F(2,3);f31=F(3,1);f32=F(3,2);f33=F(3,3);
G = reshape(vm(:, end-1), 3, 3)';
g11=G(1,1);g12=G(1,2);g13=G(1,3);g21=G(2,1);g22=G(2,2);g23=G(2,3);g31=G(3,1);g32=G(3,2);g33=G(3,3);
a = (g11*g22*g33 - g11*g23*g32 - g12*g21*g33 + g12*g23*g31 + g13*g21*g32 - g13*g22*g31);
b = (f11*g22*g33 - f11*g23*g32 - f12*g21*g33 + f12*g23*g31 + f13*g21*g32 - f13*g22*g31 - f21*g12*g33 + f21*g13*g32 + f22*g11*g33 - f22*g13*g31 ...
    - f23*g11*g32 + f23*g12*g31 + f31*g12*g23 - f31*g13*g22 - f32*g11*g23 + f32*g13*g21 + f33*g11*g22 - f33*g12*g21);
c = (f11*f22*g33 - f11*f23*g32 - f11*f32*g23 + f11*f33*g22 - f12*f21*g33 + f12*f23*g31 + f12*f31*g23 - f12*f33*g21 + f13*f21*g32 - f13*f22*g31 ...
    - f13*f31*g22 + f13*f32*g21 + f21*f32*g13 - f21*f33*g12 - f22*f31*g13 + f22*f33*g11 + f23*f31*g12 - f23*f32*g11);
d = f11*f22*f33 - f11*f23*f32 - f12*f21*f33 + f12*f23*f31 + f13*f21*f32 - f13*f22*f31;
s = roots([a b c d]);
flag = false(1,length(s));
for i = 1:length(s)
    flag(i) = isreal(s(i));
end
s = s(flag);
Fm = zeros(3,3,length(s));
for i = 1:length(s)
    Fm(:,:,i) = T2'*(F + s(i)*G)*T1;
end

end