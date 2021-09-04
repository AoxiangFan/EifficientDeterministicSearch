function [b,r] = wfit(y,x,w)
%WFIT    weighted least squares fit

% Create weighted x and y
n = size(x,2);
sw = sqrt(w);
yw = y .* sw;
xw = x .* sw(:,ones(1,n));

% Computed weighted least squares results
[b,r] = linsolve(xw,yw,struct('RECT',true));