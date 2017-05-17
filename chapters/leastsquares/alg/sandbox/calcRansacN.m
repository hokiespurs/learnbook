function N = calcRansacN(ransacconfidence,outlierprob,nPtsInSample)
% calculate the number of iterations required to have no outliers with an
% input confidence
N = log(1-ransacconfidence)./log(1-(1-outlierprob).^nPtsInSample);

end