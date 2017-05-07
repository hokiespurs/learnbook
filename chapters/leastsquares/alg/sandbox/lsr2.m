function [betacoef,R,J,CovB,MSE,ErrorModelInfo] = lsr(x,y,modelfun,varargin)
% LSR Leasr Squares Regression
% Inputs:
%   - x              : Predictor variables
%   - y              : Response values
%   - modelfun       : Model function handle @modelfun(betacoef,X)
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
%   - 4th input / 'betacoef0': Initial coefficient values
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
%   - 'ransac'               : optional structure with ransac parameters

%% Input Parsing/Checks
% Get constants from default parameters
narginchk(3,inf);
nObsEqns = numel(y);
nPredictors = numel(x);
nBetacoef = calcNbetacoef(modelfun,x); %throws error if bad modelfun

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
defaultScaleCov     = true;
defaultMaxIter      = 100;
defaultVerbose      = false;
defaultDerivstep    = eps^(1/3);
defaultRansac       = [];

% expected Values
expectedType = {'ols','wls','gls','lin','linear',...
                'nlin','nonlinear',...
                'total','tls',...
                'robust',... %still do hessian to determine linearity
                'robustlinear',...
                'robustnonlinear'};

% Input checkFunctions
checkX            = @(X) isnumeric(X) && size(X,1)==nObsEqns && all(~isnan(X(:)));
checkY            = @(X) isnumeric(X) && size(X,2)==1 && all(~isnan(X(:)));
checkModelfun     = @(X) isa(X,'function_handle');
checkBetaCoef0    = @(X) isnumeric(X) && isvector(X) && numel(X)==nBetacoef;
checkType         = @(X) any(validatestring(X, expectedType));
checkWeights      = @(X) longWeightCheck(X,nObsEqns,nPredictors); 
checkBetaCoef0Cov = @(X) isCovariance(X,nBetacoef);
checkJacobianYB   = @(X) isa(X,'function_handle');
checkJacobianYX   = @(X) isa(X,'function_handle'); 
checkRobustWgtFun = @(X) checkRobust(X,y); 
checkTune         = @(X) isnumeric(X);
checkChi2alpha    = @(X) isnumeric(X) && isscalar(X) && X>0 && X<1;
checkScaleCov     = @(X) islogical(X);
checkMaxIter      = @(X) isinteger(X) && X>0;
checkVerbose      = @(X) islogical(X);
checkDerivStep    = @(X) isnumeric(X) && isscalar(X) && X>0;
checkRansac       = @(x) isstruct(X);

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
addParameter(p,'ransac'      ,defaultRansac       ,checkRansac);

parse(p,x,y,modelfun,varargin{:});

ErrorModelInfo.parser = p; % include function input in metadata output

% get variables out of structure
betaCoef0 = p.Results.betaCoef0;
maxiter = p.Results.maxiter;
isverbose = p.Results.verbose;
scalecov = p.Results.scaleCov;
chi2alpha = p.Results.chi2alpha;
[robustTune, robustWgtFun] = getRobust(p.Results.RobustWgtFun, p.Results.Tune);
betacoef0cov = p.Results.betaCoef0Cov;
[Sx, ~] = getCovariance(p.Results.weights,nObsEqns);
if isempty(p.Results.JacobianYB)
    JybFunction = @(b,x) calcJYB(modelfun,b,x,p.Results.derivstep);
else
    JybFunction = p.Results.JacobianYB;
end
if isempty(p.Results.JacobianYX)
    JyxFunction = @(b,x) calcJYX(modelfun,b,x,p.Results.derivstep);
else
    JyxFunction = p.Results.JacobianYX;
end
%% Determine least squares type and do some checks of input values
lstype = getlstype(p.Results.type,modelfun,betaCoef0,x);

dolstypechecks(lstype,p);

[ransacparams,doransac] = checkRansac(p.Results.ransac);

if isverbose
    printPreSummary(lstype,ransacparams,p,nBetacoef);
end

%% Make function handle for Specific Least Squares Function Depending
switch lstype
    case 1
        lsrfun = @(x,y) lsrlin(x,y,modelfun,betaCoef0,Sx,betacoef0cov,...
            JybFunction,chi2alpha,scalecov,isverbose);
    case 2
        lsrfun = @(x,y) lsrnlin(x,y,modelfun,betaCoef0,Sx,betacoef0cov,...
            JybFunction,chi2alpha,scalecov,maxiter,isverbose);
    case 3
        lsrfun = @(x,y) lsrrobustlin(x,y,modelfun,betaCoef0,...
            JybFunction,robustWgtFun,robustTune,isverbose);
    case 4
        lsrfun = @(x,y) lsrrobustnlin(x,y,modelfun,betaCoef0,...
            JybFunction,robustWgtFun,robustTune,maxiter,isverbose);
    case 5
        lsrfun = @(x,y) lsrtotal(x,y,modelfun,betaCoef0,Sx,betacoef0cov,...
            JybFunction,JyxFunction,chi2alpha,scalecov,maxiter,isverbose);
end

%% Compute Least Squares (optionally use ransac)
if ~doransac
    [betacoef,R,J,CovB,MSE,ErrorModelInfo] = lsrfun(x,y);
else
    [betacoef,R,J,CovB,MSE,ErrorModelInfo] = calcransac(lsrfun,x,y,ransacparams);
end

end

function printPreSummary(lstype,p,nBetacoef)
%% Print output to the screen summarizing least squares regression params
hline = [repmat('-',1,50) '\n'];
switch lstype
    case 1
        fprintf('%sLINEAR LEAST SQUARES\n%s',hline,hline);
    case 2
        fprintf('%sNONLINEAR LEAST SQUARES\n%s',hline,hline);
    case 3
        fprintf('%sROBUST LINEAR LEAST SQUARES\n%s',hline,hline);
    case 4
        fprintf('%sROBUST NONLINEAR LEAST SQUARES\n%s',hline,hline);
    case 5
        fprintf('%sTOTAL LEAST SQUARES\n%s',hline,hline);
end

nObsEqns = numel(p.Results.y);
nPredictors = numel(p.Results.x);
ndof = nObsEqns-nBetacoef; 

fprintf('test');
end

function [betacoef,R,J,CovB,MSE,ErrorModelInfo] = lsrlin(x,y,modelfun,...
    betaCoef0,Sx,betacoef0cov,JybFunction,chi2alpha,scalecov,isverbose)
%% Calculate linear least squares
nBetacoef = calcNbetacoef(modelfun,x);
A = JybFunction(zeros(nBetacoef,1),x);
L = y;

if betacoef0asobs %add betacoef0 as observation equations
    if isverbose
       fprintf('Using betacoef0 as observation equations: yes\n'); 
    end
    
    L(end+1:end+nBetacoef)=betaCoef0;
    A = [A;eye(nBetacoef)];
    covY = blkdiag(Sx,betacoefcov);
else
    if isverbose
       fprintf('Using betacoef0 as observation equations: no\n'); 
    end    
end

W = inv(covY);

betacoef = (A'*W*A)\A'*W*L;     % unknowns
V = A * betacoef - L;           % residuals
So2 = V'*W*V/dof;               % Reference Variance
ErrorModelInfo.betacoef = betacoef;
ErrorModelInfo.V = V;
ErrorModelInfo.MSE = So2;
ErrorModelInfo.Q = inv(A'*W*A);              % cofactor
if isverbose
    fprintf('        So2        '); %
    fprintf('         betacoef(%.0f)',1:nBetacoef);
    fprintf('\n');
    fprintf('%20.10f', So2);
    fprintf('%20.10f', betacoef);
    fprintf('\n');
end
if p.Results.noscale %force it not to scale covariance
    if isverbose
        fprintf('Covariance is not scaled\n');
    end
    ErrorModelInfo.covB = inv(A'*W*A);
else
    if isverbose
        fprintf('Covariance is scaled\n');
    end
    ErrorModelInfo.covB = So2*inv(A'*W*A);
end
ErrorModelInfo.stdX = sqrt(diag(ErrorModelInfo.covB));% std of solved unknowns
ErrorModelInfo.Lhat = A * betacoef;                   % predicted L values
ErrorModelInfo.RMSE = sqrt(V'*V/nObsEqns);            % RMSE
ErrorModelInfo.r2 = var(ErrorModelInfo.Lhat)/var(L);  % R^2 Skill
% handle output variables
R = V;
MSE = So2;
CovB = ErrorModelInfo.covB;
J = A;

end

function [ransacparams,doransac] = checkRansac(ransac)
% input a possible ransac structure and return all the variables needed


end

function [betacoef,R,J,CovB,MSE,ErrorModelInfo] = calcransac(lsrfun,x,y,ransacparams)
% do ransac parameter estimation


end

function dolstypechecks(lstype,p)
%% throw a bunch of errors/warnings depending on the parameters input
% see flowchart for summary of this
%determine input variables
indsused = true(numel(p.Parameters),1);
for i=1:numel(p.UsingDefaults)
    indunused = strcmp(p.UsingDefaults{i},p.Parameters);
    indsused(indunused)=false;
end
inputParams = p.Parameters(indsused);

switch lstype
    case 1 % linear
        requiredFields = {'x'}; %error would have been thrown earlier if bad
        illegalFields = {'JacobianYX','Tune','maxiter'};
        typename = 'linear';
    case 2 % nonlinear
        requiredFields = {'betaCoef0'};
        illegalFields = {'JacobianYX','Tune'};
        typename = 'nonlinear';
    case 3 % robust linear
        requiredFields = {'RobustWgtFun'}; %error would have been thrown earlier if bad
        fprintf('robust linear\n');
        illegalFields = {'weights','betaCoef0Cov',...
                         'JacobianYX','chi2alpha','scaleCov','maxiter'};
        typename = 'robust linear';
    case 4 % robust nonlinear
        requiredFields = {'betaCoef0','RobustWgtFun'};
        illegalFields = {'weights','betaCoef0Cov','JacobianYX',...
                         'chi2alpha','scaleCov'};
        typename = 'robust nonlinear';
    case 5 % total least squares
        requiredFields = {'betaCoef0'};
        illegalFields = {'Tune'};
        typename = 'total least squares';
end
for i=1:numel(requiredFields) %% loop through all parameters NOT input
    % throw errors
    if ~any(strcmp(requiredFields{i},inputParams))
        strerr = sprintf('%s,',requiredFields{:});
        error('Missing at least one of the required parametersfor %s: %s\n',...
            typename,strerr);
    end
end
for i=1:numel(illegalFields)
    if any(strcmp(illegalFields{i},inputParams))
        strerr = sprintf('%s,',illegalFields{:});
        error('Included at least one illegal parameter for %s: %s\n',...
            typename,strerr);
    end
end
    % throw warnings
    [~, iscovariance] = getCovariance(p.Results.weights);
    if any(lstype == [1 2 5])
        if ~iscovariance && any(strcmp('chi2alpha',inputParams))
            warning('chi2alpha is meaningless without input covariance');
        elseif ~iscovariance && any(strcmp('betaCoef0Cov',inputParams))
            warning('betaCoef0Cov is meaningless without input covariance');
        end
    end
    if any(lstype == [1 3]) && ~any(strcmp('betaCoef0Cov',inputParams)) && any(strcmp('beta0',inputParams))
        warning('For linear systems, beta0 does nothing unless betaCoef0Cov is input');
    end
end

function lsrtype = getlstype(inputtype,modelfun,betaCoef0,x)
%% Determine the type of least squares
% lstype 
% 1: linear
% 2: nonlinear
% 3: robust linear
% 4: robust nonlinear
% 5: total

if isempty(inputtype) % either linear or nonlinear
    if isempty(betaCoef0) || isModelLinear(modelfun,betaCoef0,x)
        lsrtype = 1;
    else
        lsrtype = 2;
    end
elseif any(strcmp(inputtype,{'ols','wls','gls','lin','linear'}))
    lsrtype = 1;
elseif any(strcmp(inputtype,{'nlin','nonlinear'}))
    lsrtype = 2;    
elseif any(strcmp(inputtype,{'total','tls'}))
    lsrtype = 5;
elseif any(strcmp(inputtype,{'robust'})) %nonlinear and linear robust 
    if isempty(betaCoef0) || isModelLinear(modelfun,betaCoef0,x)
        lsrtype = 3; %assume no beta0 means linear
    else
        lsrtype = 4;
    end
elseif any(strcmp(inputtype,{'robustlinear'}))
    lsrtype = 3;
elseif any(strcmp(inputtype,{'robustnonlinear'}))
    lsrtype = 4;
else
    error('Unable to determine the type of least squares');
end

end

function isLinear = isModelLinear(modelfun,betacoef,x)
%% Compute the Numerical Hessian to determine linearity of modelfun
h = eps^(1/3); % optimal for central difference
Jfun = @(bn)(modelfun(bn,x(1,:)));
J = @(b) (calcPartials(Jfun,betacoef',h));
H = calcPartials(J,betacoef',h); %calculate hessian
isLinear = ~any(H(:));
end

function nBetacoef = calcNbetacoef(modelfun,x)
%% Try catch loop until the number of beta parameters doesnt throw an error
for iTestBeta=1:100 %loop and test different nBetacoef and see what works
    try
        modelfun(zeros(iTestBeta,1),x(1,:));
        nBetacoef = iTestBeta;
        break;
    catch
        % keep looping until modelfun doesnt throw an error
        % maybe I shouldnt have tried to force linear into this function
    end
end
if nBetacoef ==100
    error('Unable to determine the number of beta coefficients');
end
end

function [robustTune, robustWgtFun, isvalid] = getRobust(weightfun, tune)
%% the weightfunction can be either a string of a function
% determine which and make sure its valid
% return the tune and weightfun

robustTune = tune;
inputWgtFun = weightfun;
if ischar(inputWgtFun)
    isvalid = true;
    switch inputWgtFun
        case 'test'
            robustWgtFun = @(x)(1+x);
            defaultTune = 0.1;
        case 'test2'
            robustWgtFun = @(x)(1+x);
            defaultTune = 0.2;
        otherwise
            isvalid = false; %don't know this function
    end
    if isempty(robustTune)
        robustTune = defaultTune;
    end
else %user explicitly input a robust weight function
    robustWgtFun = inputWgtFun;
    if isempty(robustTune)
       error('User Input robustWgtFun must have a ''Tune'' Value input'); 
    end
end

end

function isvalid = checkRobust(X,y)
% if its a string, ensure its an expected value
% if its a function handle, ensure it works
if ischar(X)
    [~,~,isvalid]=getRobust(X,1); % call getRobust to see if its valid
elseif isa(X,'function_handle')
    try
        isvalid = true;
        X(y);
    catch
        isvalid = false;
    end
else
    isvalid = false;
end
end

function isvalid = longWeightCheck(weights,nObsEqns,nPredictors)
% determine if the weights input is valid
% can either be a vector of weights, or a covariance matrix
% the covariance matrix must be valid

if isvector(weights) %matrix is weights
    Sx = diag(1./weights);
else %matrix is covariance
    Sx = weights;
end

isvalid = false;

if isCovariance(Sx,nObsEqns) || isCovariance(Sx,nPredictors)
    isvalid = true;
end

end

function isvalid = isCovariance(C,ndim)
% check if a covariance matrix is valid
isvalid = false;
if nargin==2 && size(C,1)~=ndim
    %     warning('Covariance matrix must be %.0fx%.0f',ndim,ndim);
elseif ~issymmetric(C)
    %     warning('Covariance matrix must be Symmetric');
else
    [~,D]=eig(C);
    eigenvals = diag(D);
    if sum(eigenvals<0)>0
        %   warning('Covariance matrix must be positive semi-definite');
    else
        isvalid = true;
    end
end
end

function [Sx, iscovariance] = getCovariance(weights,n)
% return a covariance depending on the type of input
% vector is "weights", where covariance is inverse of those on the diagonal
% matrix is assumed to be a covariance
if isempty(weights)
    if nargin ==2
       Sx = speye(n); 
    end
    iscovariance = false;
elseif isvector(weights) %matrix is weights
    Sx = diag(1./weights);
    iscovariance = false;
else %matrix is covariance
    Sx = weights;
    iscovariance = true;
end

end

%% Numerical Partial Derivatives
function Jyb = calcJYB(modelfun,b,x,h)
%% Numerically calculate the Jacobian for modelfun(b,x) wrt b
% h (optional) is the delta for calculating the central finite difference.

if nargin==3
    h=eps^(1/3); % optimal for central difference (source: internet)
end

Jybfun = @(bn)(modelfun(bn,x));
Jyb = calcPartials(Jybfun,b',h);

end

function Jyx = calcJYX(modelfun,b,x,h)
%% Numerically calculate the Jacobian for modelfun(b,x) wrt x
% h (optional) is the delta for calculating the central finite difference.

if nargin==3
    h=eps^(1/3); % optimal for central difference (source: internet)
end

JYXfun = @(xn)(modelfun(b,xn));
Jyx = calcPartials(JYXfun,x,h);

nEqn = size(modelfun(betacoef,x),1)/size(x,1);

Jyx = bumphdiag(Jyx,nEqn);
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

function hx = bumphdiag(x,n)
    % pad a matrix with 0s while shifting n rows at a time over on blkdiag
    % This is useful when doing partial derivatives wrt predictor variables
    %
    % ie, with n=1
    %  1 2 3       1 2 3 0 0 0 0 0 0 0 0 0
    %  4 5 6  ->   0 0 0 4 5 6 0 0 0 0 0 0
    %  7 8 9       0 0 0 0 0 0 7 8 9 0 0 0
    %  2 4 6       0 0 0 0 0 0 0 0 0 2 4 6
    %
    % ie, with n=2
    %  1 2 3       1 2 3 0 0 0
    %  4 5 6  ->   4 5 6 0 0 0
    %  7 8 9       0 0 0 7 8 9
    %  2 4 6       0 0 0 2 4 6
    %
    [M,N]=size(x);
    irow = kron((1:M)',ones(1,N));
    icol = nan(N,M/n);
    icol(:) = 1:numel(icol);
    icol = kron(icol.',ones(n,1));
    hx = sparse(irow,icol,x);
end