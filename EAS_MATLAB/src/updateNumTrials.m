function maxTrials = updateNumTrials(oneOverNPts, logOneMinusConf, numCurInliers, maxTrials,k)

ratioOfInliers = numCurInliers * oneOverNPts;
if ratioOfInliers > 1 - eps('double')
   newNum = 1;
else
  ratio_t = ratioOfInliers^k;
  if ratio_t > eps('double')
    logOneMinusRatio_t = log(1 - ratio_t);
    newNum = ceil(logOneMinusConf / logOneMinusRatio_t);
  else
    newNum = intmax('int32');
  end
end

if maxTrials > newNum
   maxTrials = newNum;
end
end