function d = SampsonDistanceF(pts1h, pts2h, f)
pfp = (pts2h' * f)';
pfp = pfp .* pts1h;
d = sum(pfp, 1) .^ 2;
epl1 = f * pts1h;
epl2 = f' * pts2h;
d = d ./ (epl1(1,:).^2 + epl1(2,:).^2 + epl2(1,:).^2 + epl2(2,:).^2);
d = sqrt(d);
end