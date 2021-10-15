function [Dx, Dy] = preSampsonDistanceH_all(X1,X2)

N = size(X1,2);
Dx = [X1'.* repmat(X2(3,:)',1,3), zeros(N,3), -X1'.* repmat(X2(1,:)',1,3)];
Dy = [zeros(N,3), X1'.* repmat(X2(3,:)',1,3) , -X1'.* repmat(X2(2,:)',1,3)];
