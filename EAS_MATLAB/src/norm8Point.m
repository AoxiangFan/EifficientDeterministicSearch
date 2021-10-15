function f = norm8Point(X, Y)
% Normalize the points
num = size(X, 2);

[X, T1] = normalise2dpts(X);
[Y, T2] = normalise2dpts(Y);
if any(any(isnan(X))) || any(any(isnan(Y)))
    f = rand(3,3);
    return
end

% Compute the constraint matrix
m = zeros(num, 9);
for idx = 1: num
  m(idx,:) = [...
    X(1,idx)*Y(1,idx), X(2,idx)*Y(1,idx), Y(1,idx), ...
    X(1,idx)*Y(2,idx), X(2,idx)*Y(2,idx), Y(2,idx), ...
                 X(1,idx),              X(2,idx), 1];
end

% Find out the eigen-vector corresponding to the smallest eigen-value.
[~, ~, vm] = svd(m, 0);
f = reshape(vm(:, end), 3, 3)';

% Enforce rank-2 constraint
[u, s, v] = svd(f);
s(end) = 0;
f = u * s * v';

f = f/norm(f(:));

f = T2'*f*T1;

end