function D = GenerateEmbeddings(X, Y, option)

X0 = X;
Y0 = Y;
N = size(X,1);
 
[Xn,~,~] = DataNorm(X0);
[Yn,~,~] = DataNorm(Y0);
Xt = [Xn';ones(1,N)];
Yt = [Yn';ones(1,N)];

if strcmp(option, 'F')
    D = zeros(9,N);
    for nn = 1:N
        D(:,nn) = kron(Yt(:,nn),Xt(:,nn));
    end
    nml = sqrt(sum(D.^2));
    D = D*diag(1./nml);
elseif strcmp(option, 'H')
    D = [];
    ooo  = zeros(1,3);
    for k=1:N
        p1 = Xt(:,k);
        p2 = Yt(:,k);
        D = [ D;
            p1'*p2(3) ooo -p1'*p2(1)
            ooo p1'*p2(3) -p1'*p2(2)
            ];
    end
    D = D';
    nml = sqrt(sum(D.^2));
    D = D*diag(1./nml);
elseif strcmp(option, 'At')
    D = [];
    ooo  = zeros(1,3);
    for k=1:N
        p1 = Xt(:,k);
        p2 = Yt(:,k);
        D = [ D;
            p1' ooo -p2(1)
            ooo p1' -p2(2)
            ];
    end
    D = D';
    nml = sqrt(sum(D.^2));
    D = D*diag(1./nml);
elseif strcmp(option, 'A')
    D = [Xn,Yn,ones(N,1)]';
    nml = sqrt(sum(D.^2));
    D = D*diag(1./nml);
end