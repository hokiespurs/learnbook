function indSubset = getRansacSubset(nTotal,nDesired,nElements,randseed)
% calculate random subset of the data as a nDesired x nElements index array
if nargin==4
   rng(randseed); 
end
nTotalPossible = nchoosek(nTotal,nElements);
if nDesired>nTotalPossible
   error('impossible to get that number of unique elements'); 
end

indSubset = nan(nDesired,nElements);
iRow = 0;

while iRow<nDesired
    testSubset = sort(randperm(nTotal,nElements));
    if ~any(all(testSubset==indSubset,2))
        iRow = iRow + 1;
        indSubset(iRow,:)=testSubset;
    end
end

end

%%

function [ransacparams,doransac] = getRansac(ransac,nBetacoef,nObservations,nObsEqnPerObservation)
% input a possible ransac structure and return all the variables needed
% Need to add more robust checks here
if isempty(ransac)
    doransac = false;
    ransacparams = [];
else
    doransac = true;
    if ~isfield(ransac,'confidence')
        ransac.confidence = 0.95;
    end
    if ~isfield(ransac,'nSamplePoints')
        ransac.nSamplePoints = nBetacoef;
    else
        if ransac.nSamplePoints < nBetacoef
            error('The number of points per sample is smaller than the number of beta coefficients');
        end
    end
    if ~isfield(ransac,'outlierpercent')
        ransac.outlierpercent = 0.10; % 10 percent outliers
    end
    if floor(nObservations*(1-ransac.outierpercent)) * nObsEqnPerObservation < nBetacoef
        error('expected outliers too large, number of inliers is an underconstrained solution');
    end
    if ~isfield(ransac,'residualThresh')
        ransac.residualThresh = 1;
    end
    if ~isfield(ransac,'iterations')
        ransac.iterations = 1;
    end
    if ~isfield(ransac,'earlytermination')
        ransac.iterations = false;
    end
    
    ransacparams = ransac;
    ransacparams.unroundedN = calcRansacN(ransacparams.confidence,ransacparams.outlierpercent,nObservations);
    ransacparams.N = ceil(ransacparams.unroundedN);
end
   
end
