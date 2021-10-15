function [bestF, bestInliers] = H2F(X, Y, H, th)

MAX_SAM = 100;
sam_sizH = 6;
sam_sizO = 4;

d = SampsonDistanceH_all(X, Y, H);
hinl = find(d < 4*th);
nhinl = find(d > 5*th);

if (length(hinl) < sam_sizH) || (length(nhinl) < sam_sizO)
    bestF = ones(3,3);
    bestInliers = [];
    return;
end

Xo = X(:, nhinl);
Yo = Y(:, nhinl);

Xs = H*Xo;
Ys = Yo;

Xh = X(:, hinl);
Yh = Y(:, hinl);

len_xs = size(Xs, 2);
m_i = sam_sizO;
max_sam = MAX_SAM;
no_sam = 0;

numBestInliers = 0;
bestInliers = [];
bestF = [];

while no_sam <= max_sam
    sample = randperm(len_xs, 2);
    no_sam = no_sam + 1;
    ec = cross(cross(Xs(:, sample(1)), Ys(:, sample(1))), cross(Xs(:, sample(2)), Ys(:, sample(2))));
    ec = ec/ec(3);
    aFt = skew_sym(ec/norm(ec));
    Ftmp = aFt*H;
    d = SampsonDistanceF(Xo, Yo, Ftmp);
    inliersO = find(d<=2*th);
    no_i = length(inliersO);
    if no_i > m_i
        m_i = no_i;
        max_sam = updateMaxTrials(no_i, max_sam, length(nhinl), 0.999999, 2);
        [F, inliers] = innerFH(Xh, Yh, Xo(:, inliersO), Yo(:, inliersO), X, Y, th, 15, sam_sizH, sam_sizO);
        if length(inliers) > numBestInliers
            numBestInliers = length(inliers);
            bestInliers = inliers;
            bestF = F;
        end
    end
end
    








    