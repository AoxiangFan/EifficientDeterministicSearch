function [B, distances] = DPCP_solver(Xtilde,c,mu_min,maxiter)

% This code slightly changed from release of NeuRIPS2018 paper "Dual principal component pursuit: Improved analysis and efficient algorithms".
if nargin < 3
    mu_min = 1e-15; maxiter = 200;
end

if nargin < 4
    maxiter = 200;
end

mu_0 = 1e-2; alpha = 1e-3; beta = 1/2;

[D, N] = size(Xtilde);
obj = @(b)norm(Xtilde'*b,1);

% initialization
% [B_0,~] = eigs(Xtilde*Xtilde',c,'SM');
[B_0, diag_0] = eig(Xtilde*Xtilde');
[~, ind] = sort(diag(diag_0));
B_0 = B_0(:, ind(1:c));
%bo = normc(randn(D,1));

for j = 1:c
    i = 1;
    b = B_0(:,j);
    mu = mu_0;
    if j == 1
        obj_old = obj(b);
        while mu>mu_min & i<= maxiter
            i = i+1;
            grad = sum( sign(b'*Xtilde).*Xtilde, 2);
            grad_norm = norm(grad)^2;
            %%% line search
            bk = b - mu*grad;
            while (obj( bk ./ norm(bk) )> obj_old - alpha*mu*grad_norm)& mu>mu_min;
                mu = mu*beta;
                bk = b - mu*grad;
            end
            b = bk ./ norm(bk);obj_old = obj(b);
        end
    else
        b = normc(b - B*(B'*b)); obj_old = obj(b);
        while mu>mu_min && i<= maxiter
            i = i+1;
            
            %%% line search
            grad = sum( sign(b'*Xtilde).*Xtilde, 2);
            % grad = grad - B*(B'*b);
            grad_norm = norm(grad)^2;
            bk = b - mu*grad;
            while (obj( bk ./ norm(bk) )> obj_old - alpha*mu*grad_norm) && mu>mu_min
                mu = mu*beta;
                bk = b - mu*grad;
            end
            bk = bk - B*(B'*bk);
            b = bk ./ norm(bk);
            obj_old = obj(b);
        end
    end
    B(:,j) = b;
    
end

distances = zeros(1,N);
for j = 1 : N
    distances(j) = norm(B' * Xtilde(:,j));
end

end
