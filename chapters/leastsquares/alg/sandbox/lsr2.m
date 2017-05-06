function [betacoef,R,J,CovB,MSE,ErrorModelInfo] = lsr(x,y,modelfun,varargin)
% LSR 
% Inputs:
%   - x              : Predictor variables
%   - y              : Response values
%   - modelfun       : Model function handle @modelfun(betacoef,X)
%   - betacoef0      : Initial coefficient values
%   - options        : Structure with optional parameters
% 
% Outputs:
%   - betacoef       : Estimated regression coefficients
%   - R              : Residuals
%   - J              : Jacobian
%   - CovB           : Estimated Variance Covariance Matrix
%   - MSE            : Mean Squared Error (Computed Reference Variance)
%   - ErrorModelInfo : Information about error model
% 
% Optional Parameters:
%   - 'type'                 : 'linear/ols/wls/gls','nonlinear'(default),'tls'
%   - 'beta0'                : required vector for nonlinear and tls
%   - 'stochastic'           : [empty] (default), vector (weights), covariance matrix (S)
%   - 'AnalyticalJacobian'   : [empty] (default), optional Jfun(beta,x)
%   - 'AnalyticalBfunction'  : [empty] (default), optional Bfun(beta,x)
%   - 'noscale'              : scale covariance...false (default), true 
%   - 'beta0Covariance'      : [empty] (default), (m x m) covariance
%   - 'betacoef0asObs'       : false (default), option to use guess in modelfun
%   - 'chi2alpha'            : 0.05 (default), alpha value for confidence
%   - 'RobustWgtFun'         : Robust Weight Function
%   - 'Tune'                 : Robust Wgt Tuning Function

%% Input Parsing/Checks
% Get constants from default parameters
narginchk(3,inf);
nObsEqns = numel(y);
nPredictors = numel(x);
nBetacoef = calcNbetacoef(modelfun,x); %throws error if bad modelfun
ndof = nObsEqns-nBetacoef; 

% default values
defaultBetaCoef0    = [];
defaultType         = [];
defaultWeights      = [];
defaultBetaCoef0Cov = [];
defaultJacobianYB   = [];
defaultJacobianYX   = [];
defaultRobustWgtFun = [];
defaultTune         = false;
defaultchi2alpha    = 0.05;
defaultScaleCov   = [];
defaultMaxIter      = 100;
defaultVerbose      = false;
defaultDerivstep    = eps^(1/3);

% expected Values
expectedType = {'ols','wls','gls','lin','linear',...
                'nlin','nonlinear',...
                'total','tls',...
                'robust',... %still do hessian to determine linearity
                'robustlinear',...
                'robustnonlinear'};

% Input checkFunctions
checkX            = @(X) isnumeric(X) && size(X,1)==nObsEqns;
checkY            = @(X) isnumeric(X) && size(X,2)==1;
checkModelfun     = @(X) isa(X,'function_handle');
checkBetaCoef0    = @(X) isnumeric(X) && isvector(X) && numel(X)==nObsEqns;
checkType         = @(X) any(validatestring(X, expectedType));
checkWeights      = @(X) longWeightCheck(X,nObsEqns,nPredictors); 
checkBetaCoef0Cov = @(X) longBetaCoef0covCheck(X,nBetacoef);
checkJacobianYB   = @(X) isa(X,'function_handle');
checkJacobianYX   = @(X) isa(X,'function_handle'); 
checkRobustWgtFun = @(X) checkRobust(xparam,y); 
checkTune         = @(X) isnumeric(X);
checkChi2alpha    = @(X) isnumeric(X) && isscalar(X) && X>0 && X<1;
checkScaleCov     = @(X) islogical(X);
checkMaxIter      = @(X) isinteger(X) && X>0;
checkVerbose      = @(X) islogical(X);
checkDerivStep    = @(X) isnumeric(X) && isscalar(X) && X>0;

% parse inputs
p = inputParser;

addRequired(p,'x'                                 ,checkX);
addRequired(p,'y'                                 ,checkY);
addRequired(p,'modelfun'                          ,checkModelfun);

addOptional(p,'betaCoef0'    ,defaultBetaCoef0    ,checkBetaCoef0);

addParameter(p,'type'        ,defaultType         ,checkType);
addParameter(p,'weights'     ,defaultWeights      ,checkWeights);
addParameter(p,'betaCoef0Cov',defaultBetaCoef0Cov ,checkBetaCoef0Cov);
addParameter(p,'JacobianYB'  ,defaultJacobianYB   ,checkJacobianYB);
addParameter(p,'JacobianYX'  ,defaultJacobianYX   ,checkJacobianYX);
addParameter(p,'RobustWgtFun',defaultRobustWgtFun ,checkRobustWgtFun);
addParameter(p,'Tune'        ,defaultTune         ,checkTune);
addParameter(p,'chi2alpha'   ,defaultchi2alpha    ,checkChi2alpha);
addParameter(p,'scaleCov'    ,defaultScaleCov     ,checkScaleCov);
addParameter(p,'maxiter'     ,defaultMaxIter      ,checkMaxIter);
addParameter(p,'verbose'     ,defaultVerbose      ,checkVerbose);
addParameter(p,'derivstep'   ,defaultDerivstep    ,checkDerivStep);

parse(p,x,y,modelfun,varargin{:});

ErrorModelInfo.parser = p; % include function input in metadata output

% get variables out of structure
maxiter = p.Results.maxiter;
isverbose = p.Results.verbose;
scalecov = p.Results.scaleCov;
chi2alpha = p.Results.chi2alpha;
robustTune = p.Results.Tune;
robustWgtFun = p.results.RobustWgtFun;
betacoef0cov = p.Results.betaCoef0Cov;
[Sx, iscovariance, isweights] = getCovariance(p.Results.weights);

JybFunction = @(b,x) calcJYB(modelfun,b,x,p.Results.derivstep);
JyxFunction = @(b,x) calcJYX(modelfun,b,x,p.Results.derivstep);


%% Call Specific Least Squares Function Depending on Computed Type
% lstype 
% 1: linear
% 2: nonlinear
% 3: robust linear
% 4: robust nonlinear
% 5: total
lstype = getlstype(p);


end

function isvalid = longWeightCheck(X,nObsEqns,nPredictors)
% if weights
% no negatives
% must be the right length
% 
% if covariance
% no negative on diagonal
% positive semi definite
% symmetric
% must be the right size

isvalid = true;
end

function isvalid = longBetaCoef0covCheck(X,nBetacoef)

isvalid = true;
end

function [Sx, iscovariance, isweights] = getCovariance(weights)

if isvector(weights)
    Sx = diag(weights);
    iscovariance = false;
    isweights = true;
elseif ismatrix
    Sx = weights;
    iscovariance = true;
    isweights = true;
end

end

function dfdxn = calcPartials(f,x,h)
% f      : model as a function of one vector @f(x)
% xi     : values of x to evaluate partial at
% n      : which partial to calculate
% h      : what step increment for finite differencing
%
% dfdxn  : partial f wrt x_n
nObservations = numel(f(x));
nVariables = size(x,2);
dfdxn= nan(nObservations,nVariables);
for i=1:nVariables
    %central derivative
    ix1 = x;
    ix1(:,i)=ix1(:,i)-h/2;
    
    ix2 = x;
    ix2(:,i)=ix2(:,i)+h/2;
    %evaluate slope
    dfdxn(:,i) = (f(ix2)-f(ix1))/h; 
end

end

function Jyb = calcJYB(modelfun,b,x,h)

if nargin==3
    h=eps^(1/3); %optimal for central difference
end
Jybfun = @(bn)(modelfun(bn,x));
Jyb = calcPartials(Jybfun,b',h);

end

function Jyx = calcJYX(modelfun,b,x,h)

if nargin==3
    h=eps^(1/3); %optimal for central difference
end

JYXfun = @(xn)(modelfun(b,xn));
Jyx = calcPartials(JYXfun,x,h);

nEqn = size(modelfun(betacoef,x),1)/size(x,1);

Jyx = bumphdiag(Jyx,nEqn);
end

function hx = bumphdiag(x,n)
    
    [M,N]=size(x);
    irow = kron((1:M)',ones(1,N));
    icol = nan(N,M/n);
    icol(:) = 1:numel(icol);
    icol = kron(icol.',ones(n,1));
    hx = sparse(irow,icol,x);
end