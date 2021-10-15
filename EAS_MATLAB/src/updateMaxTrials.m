function maxTrials = updateMaxTrials(numCurInliers, maxTrials, N, confidence, k)

ratioOfInliers = numCurInliers/N;
if ratioOfInliers > 1 - 1e-16
    newNum = 1;
else
    ratio_t = ratioOfInliers^k;
    if ratio_t > 1e-16
        logOneMinusRatio_t = log(1 - ratio_t);
        newNum = ceil(log(1 - confidence)/logOneMinusRatio_t);
    else
        newNum = 1e16;
    end
end
    
if maxTrials > newNum
    maxTrials = newNum;
end

end