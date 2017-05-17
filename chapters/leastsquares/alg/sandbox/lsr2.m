function [betacoef,R,J,CovB,MSE,ErrorModelInfo] = lsr2(x,y,modelfun,varargin)
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
%   - 'weights'              : [empty] (default), vector (weights), covariance matrix (S)
%   - 'AnalyticalJacobian'   : [empty] (default), optional Jfun(beta,x)
%   - 'AnalyticalBfunction'  : [empty] (default), optional Bfun(beta,x)
%   - 'noscale'              : scale covariance...false (default), true 
%   - 'beta0Covariance'      : [empty] (default), (m x m) covariance
%   - 'betacoef0asObs'       : false (default), option to use guess in modelfun
%   - 'chi2alpha'            : 0.01 (default), alpha value for confidence
%   - 'RobustWgtFun'         : Robust Weight Function
%   - 'Tune'                 : Robust Wgt Tuning Function

%% Input Parsing/Checks
% Get constants from default parameters
narginchk(3,inf);
nObsEqns = numel(y);
nObservations = size(x,1);
nObsEqnPerObservation = nObsEqns/nObservations;
nPredictors = numel(x);
nBetacoef = calcNbetacoef(modelfun,x); %throws error if bad modelfun
ndof = nObsEqns - nBetacoef;

% default values
defaultBetaCoef0    = [];
defaultType         = [];
defaultWeights      = [];
defaultBetaCoef0Cov = [];
defaultJacobianYB   = [];
defaultJacobianYX   = [];
defaultRobustWgtFun = [];
defaultTune         = [];
defaultchi2alpha    = 0.05;
defaultScaleCov     = true;
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
checkMaxIter      = @(X) mod(X,1) == 0 && X>0;
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
betaCoef0 = p.Results.betaCoef0(:);
maxiter = p.Results.maxiter;
isverbose = p.Results.verbose;
scalecov = p.Results.scaleCov;
chi2alpha = p.Results.chi2alpha;
[robustTune, robustWgtFun] = getRobust(p.Results.RobustWgtFun, p.Results.Tune);
betacoef0cov = p.Results.betaCoef0Cov;
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
lstype = getlstype(p.Results.type,modelfun,betaCoef0,x,p);

dolstypechecks(lstype,p);

[Sx, isinputcovariance] = getCovariance(lstype,p.Results.weights,x);

if isverbose
    printPreSummary(lstype,p,nBetacoef,nObsEqns);
end

%% Make function handle for Specific Least Squares Function Depending
switch lstype
    case 1
        lsrfun = @(x,y) lsrlin(x,y,modelfun,betaCoef0,Sx,betacoef0cov,...
            JybFunction,scalecov,isverbose);
    case 2
        lsrfun = @(x,y) lsrnlin(x,y,modelfun,betaCoef0,Sx,betacoef0cov,...
            JybFunction,[],scalecov,maxiter,isverbose);
    case 3
        lsrfun = @(x,y) lsrrobustlin(x,y,modelfun,betaCoef0,betacoef0cov,...
            JybFunction,scalecov,isverbose,robustWgtFun,robustTune);
    case 4
        lsrfun = @(x,y) lsrrobustnlin(x,y,modelfun,betaCoef0,...
            JybFunction,robustWgtFun,robustTune,maxiter,isverbose);
    case 5
        lsrfun = @(x,y) lsrnlin(x,y,modelfun,betaCoef0,Sx,betacoef0cov,...
            JybFunction,JyxFunction,scalecov,maxiter,isverbose);
end

%% Compute Least Squares
[betacoef,R,J,CovB,MSE,ErrorModelInfo] = lsrfun(x,y);

%% Do Chi2 Test if covariance and linear/nonlinear
if isinputcovariance  
    chi2 = ndof * MSE;
    chi2low  = chi2inv(chi2alpha/2  , ndof);
    chi2high = chi2inv(1-chi2alpha/2, ndof);
    if chi2low<chi2 && chi2<chi2high
       ErrorModelInfo.chi2.pass = true;
       passfailstr = 'PASS';
    else
        ErrorModelInfo.chi2.pass = false;
        passfailstr = 'FAIL';
    end
    ErrorModelInfo.chi2.alpha          = chi2alpha;
    ErrorModelInfo.chi2.calculatedchi2 = chi2;
    ErrorModelInfo.chi2.chi2low        = chi2low;
    ErrorModelInfo.chi2.chi2high       = chi2high;
    if isverbose
        fprintf('\n\tCHI^2 GOODNESS OF FIT TEST (Significance=%.2f, dof=%.0f, So2=%.3f)\n',chi2alpha,ndof,MSE);
        fprintf('\t  Ho %.3f == 1\n',MSE);
        fprintf('\t  H1 %.3f =/= 1\n',MSE);
        fprintf('\t  Statistical Test : (%.3f < %.3f < %.3f)\n',chi2low/ndof,chi2/ndof,chi2high/ndof);
        if scalecov && ErrorModelInfo.chi2.pass % scale but it passed test 
            fprintf('\t  **%s**  * Test Passed, consider setting ''scalecov'' == false\n',passfailstr);
        elseif ~scalecov && ~ErrorModelInfo.chi2.pass % no scale but didnt pass
            fprintf('\t  **%s**  * Test Failed, consider setting ''scalecov'' == true\n',passfailstr);
        else
            fprintf('\t  **%s**  \n',passfailstr);
        end
    end
end

end

function printPreSummary(lstype,p,nBetacoef,nObsEqns)
%% Print output to the screen summarizing least squares regression params
nObsEqns = numel(p.Results.y);
nPredictors = numel(p.Results.x);
ndof = nObsEqns-nBetacoef; 

hline = repmat('-',1,50);
switch lstype
    case 1
        fprintf('%s\nLINEAR LEAST SQUARES\n%s\n',hline,hline);
    case 2
        fprintf('%s\nNONLINEAR LEAST SQUARES\n%s\n',hline,hline);
    case 3
        fprintf('%s\nROBUST LINEAR LEAST SQUARES\n%s\n',hline,hline);
    case 4
        fprintf('%s\nROBUST NONLINEAR LEAST SQUARES\n%s\n',hline,hline);
    case 5
        fprintf('%s\nTOTAL LEAST SQUARES\n%s\n',hline,hline);
end

fprintf('\t # of Observation Equations : %.0f\n',nObsEqns);
fprintf('\t # of Predictor Variables   : [%.0f x %.0f]\n',size(p.Results.x));
fprintf('\t # of Beta Coefficients     : %.0f\n',nBetacoef);
fprintf('\t # of Degrees of Freedom    : %.0f\n',ndof);

end

function [betacoef,R,J,CovB,MSE,ErrorModelInfo] = lsrrobustlin(x,y,...
    modelfun,betaCoef0,betacoef0cov,JybFunction,scalecov,isverbose,...
    robustWgtFun,robustTune)
%% Calculate Robust Linear Least Squares
if isverbose
    n = calcNbetacoef(modelfun,x);  % number of unknowns
    fprintf('LINEAR ROBUST\n');
    fprintf('\niter :        So2        ');
    fprintf('         betacoef(%.0f)',1:n);
    fprintf('\n');
end
%initialize while loop
dMSE = 1;
lastMSE = inf;
iter = 0;
MSETHRESH = 1e-10;
MAXROBUSTITER = 200;

W = ones(numel(y),1);
while iter<MAXROBUSTITER && abs(dMSE)>MSETHRESH
    Sx = diag(1./W);
    [betacoef,R,J,CovB,MSE,ErrorModelInfo] = lsrlin(x,y,modelfun,...
    betaCoef0,Sx,betacoef0cov,JybFunction,scalecov,false);
    dMSE = lastMSE-MSE;
    lastMSE = MSE;
    iter = iter+1;
%     SigmaR = std(R); % Not Robust to Outliers
    SigmaR = median(abs(R))*0.6745; %estimate std using median
    %leverage and hat in robust least squares according to matlab
    % Weights points low that are outliers in the x dimension
    % https://www.mathworks.com/help/stats/robustfit.html
    % https://www.mathworks.com/help/stats/hat-matrix-and-leverage.html
    H = x*inv(x.'*x)*x.'; %Hat Matrix
    h = diag(H); %leverage
    
    Rnormalized = R./(robustTune * SigmaR .* sqrt(1-h));
    W = robustWgtFun(Rnormalized);
    W(W<=eps)=eps;
    if isverbose
        fprintf('%3.0f : ',iter);
        fprintf('%20.10f', MSE);
        fprintf('%20.10f', betacoef);
        fprintf('\n');
    end
end

end

function [betacoef,R,J,CovB,MSE,ErrorModelInfo] = lsrrobustnlin(x,y,modelfun,betaCoef0,...
            JybFunction,robustWgtFun,robustTune,maxiter,isverbose)
%% Calculate Robust Nonlinear Least Squares
scalecov = false;
JyxFunction = [];
betacoef0cov = [];

if isverbose
    n = calcNbetacoef(modelfun,x);  % number of unknowns
    fprintf('LINEAR ROBUST\n');
    fprintf('\niter :        So2        ');
    fprintf('         betacoef(%.0f)',1:n);
    fprintf('\n');
end
%initialize while loop
dMSE = 1;
lastMSE = inf;
iter = 0;
MSETHRESH = 1e-10;
MAXROBUSTITER = 200;
W = ones(numel(y),1);
while iter<MAXROBUSTITER && abs(dMSE)>MSETHRESH
    Sx = diag(1./W);
    [betacoef,R,J,CovB,MSE,ErrorModelInfo] = lsrnlin(x,y,modelfun,betaCoef0,Sx,betacoef0cov,...
            JybFunction,JyxFunction,scalecov,maxiter,false);
    dMSE = lastMSE-MSE;
    lastMSE = MSE;
    iter = iter+1;
%     SigmaR = std(R); % Not Robust to Outliers
    SigmaR = median(abs(R))*0.6745; %estimate std
    H = x*inv(x.'*x)*x.'; %Hat Matrix
    h = diag(H); %leverage
    %leverage and hat in robust least squares according to matlab
    % https://www.mathworks.com/help/stats/robustfit.html
    % https://www.mathworks.com/help/stats/hat-matrix-and-leverage.html
    Rnormalized = R./(robustTune * SigmaR .* sqrt(1-h));
    W = robustWgtFun(Rnormalized);
    W(W<=eps)=eps;
    if isverbose
        fprintf('%3.0f : ',iter);
        fprintf('%20.10f', MSE);
        fprintf('%20.10f', betacoef);
        fprintf('\n');
    end
end

end

function [betacoef,R,J,CovB,MSE,ErrorModelInfo] = lsrnlin(x,y,modelfun,betaCoef0,Sx,betacoef0cov,...
            JybFunction,JyxFunction,scalecov,maxiter,isverbose)
%% Calculate Nonlinear and Total Least Squares
    m = numel(y);                   % number of observations
    n = calcNbetacoef(modelfun,x);  % number of unknowns
    dof = m-n;                      % degrees of freedom

    if isempty(JyxFunction)
        istls = false;
    else
        istls = true;
    end

    if isverbose
        fprintf('\niter :        So2        ');
        fprintf('         betacoef(%.0f)',1:n);
        fprintf('\n');
    end

    %initialize while loop parameters
    So2 = inf;
    dSo2 = 1;
    iter = 0;
    betacoef = betaCoef0;                  % set first guess at unknowns
    while dSo2>0 && iter<maxiter %loop until So2 increases or exceed maxiter iterations
            J = JybFunction(betacoef,x);
            K = calcK(modelfun,betacoef,x,y);
            if istls
                B = JyxFunction(betacoef,x);
                W = inv(B*Sx*B');            %equivalent weight matrix
            else
                W = inv(Sx);
            end
            
            if ~isempty(betacoef0cov)
                if istls
                    [J,K,W,B]=addBetacoef0TLS(J,K,B,Sx,betaCoef0,betacoef0cov,betacoef);
                else
                    [J,K,W]=addBetacoef0Nlin(J,K,Sx,betaCoef0,betacoef0cov,betacoef);
                end
            end
        
        dbetacoef = (J'*W*J)\J'*W*K;       % Loop Delta Estimate
        betacoef = betacoef + dbetacoef;   % Loop Estimate
        V = K;                             % Residuals
        dSo2 = So2 - V'*W*V/dof;           % Change in Reference Variance
        So2 = (V'*W*V)/dof;                % Reference Variance
        iter = iter + 1;
        % print status to screen
        if isverbose
            fprintf('%3.0f : ',iter);
            fprintf('%20.10f', So2);
            fprintf('%20.10f', betacoef);
            fprintf('\n');
        end        
    end
    if istls
        if ~isempty(betacoef0cov)
            covX = blkdiag(Sx,betacoefcov);
        else
            covX = Sx;
        end
       Robs = covX * B' * W * V;
       ErrorModelInfo.Robs = Robs;
    end
    Q = inv(J'*W*J);              % cofactor
    if scalecov                   % option to not scale covariance
       Sx = So2 * Q;
    else
       Sx = Q;
    end
    stdX = sqrt(diag(Sx));        % std of solved unknowns
    Lhat = J * betacoef;          % predicted L values
    RMSE = sqrt(V'*V/m);          % RMSE
    % handle output variables
    R = V;
    CovB = Sx;
    MSE = So2;

    ErrorModelInfo.m = m;
    ErrorModelInfo.n = n;
    ErrorModelInfo.dof = dof;
    ErrorModelInfo.betacoef = betacoef;
    ErrorModelInfo.R = V;
    ErrorModelInfo.So2 = So2;
    ErrorModelInfo.Q = Q;
    ErrorModelInfo.CovB = Sx;
    ErrorModelInfo.stdX = stdX;
    ErrorModelInfo.Lhat = Lhat;
    ErrorModelInfo.r2 = NaN; %no r2 for nonlinear
    ErrorModelInfo.RMSE = RMSE;        
        
end

function [J,K,W,B]=addBetacoef0TLS(J,K,B,covX,betacoef0,betacoefcov,betacoef)
%adds beta coefficient 
Jadd = eye(numel(betacoef0));
J = [J; Jadd];

K = [K; betacoef0-betacoef];

B = blkdiag(B,-eye(numel(betacoef0)));

covX = blkdiag(covX,betacoefcov);

W = inv(B*covX*B');

end

function [J,K,W]=addBetacoef0Nlin(J,K,covY,betacoef0,betacoefcov,betacoef)
Jadd = eye(numel(betacoef0));
J = [J; Jadd];

K = [K; betacoef0-betacoef];

covY = blkdiag(covY,betacoefcov);

W = inv(covY);

end

function [betacoef,R,J,CovB,MSE,ErrorModelInfo] = lsrlin(x,y,modelfun,...
    betaCoef0,Sx,betacoef0cov,JybFunction,scalecov,isverbose)
%% Calculate linear least squares
nBetacoef = calcNbetacoef(modelfun,x);
A = JybFunction(zeros(nBetacoef,1),x);
L = y;

if ~isempty(betacoef0cov) %add betacoef0 as observation equations
    L(end+1:end+nBetacoef)=betaCoef0;
    A = [A;eye(nBetacoef)];
    Sx = blkdiag(Sx,betacoef0cov);
end

%do least squares calculation
W = inv(Sx);
m = numel(L);                   % number of observations
n = nBetacoef;                  % number of unknowns
dof = m-n;                      % degrees of freedom
X = (A'*W*A)\A'*W*L;            % unknowns
V = A * X - L;                  % residuals
So2 = V'*W*V/dof;               % Reference Variance
Q = inv(A'*W*A);                % cofactor
if scalecov                     % option to not scale covariance
   Sx = So2 * Q;
else
   Sx = Q;
end
stdX = sqrt(diag(Sx));          % std of solved unknowns
Lhat = A * X;                   % predicted L values
r2 = var(Lhat)/var(L);          % R^2 Skill
RMSE = sqrt(V'*V/m);            % RMSE

%assemble output variables and structure
betacoef = X;
R = V;
J = A;
CovB = Sx;
MSE = So2;

ErrorModelInfo.m=m;
ErrorModelInfo.n=n;
ErrorModelInfo.dof=dof;
ErrorModelInfo.betacoef=X;
ErrorModelInfo.V=V;
ErrorModelInfo.So2=So2;
ErrorModelInfo.Q=Q;
ErrorModelInfo.CovB=Sx;
ErrorModelInfo.stdX=stdX;
ErrorModelInfo.Lhat=Lhat;
ErrorModelInfo.r2=r2;
ErrorModelInfo.RMSE=RMSE;         

% print info to screen if verbose
if isverbose
    fprintf('        So2        '); %
    fprintf('         betacoef(%.0f)',1:nBetacoef);
    fprintf('\n');
    fprintf('%20.10f', So2);
    fprintf('%20.10f', betacoef);
    fprintf('\n');
end
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
    [~, isinputcovariance] = getCovariance([],p.Results.weights,[]);
    if any(lstype == [1 2 5])
        if ~isinputcovariance && any(strcmp('chi2alpha',inputParams))
            warning('chi2alpha is meaningless without input covariance');
        elseif ~isinputcovariance && any(strcmp('betaCoef0Cov',inputParams))
            warning('betaCoef0Cov is meaningless without input covariance');
        end
    end
    if any(lstype == [1 3]) && ~any(strcmp('betaCoef0Cov',inputParams)) && any(strcmp('beta0',inputParams))
        warning('For linear systems, beta0 does nothing unless betaCoef0Cov is input');
    end
end

function lsrtype = getlstype(inputtype,modelfun,betaCoef0,x,p)
%% Determine the type of least squares
% lstype 
% 1: linear
% 2: nonlinear
% 3: robust linear
% 4: robust nonlinear
% 5: total
isRobust = ~any(strcmp(p.UsingDefaults,'RobustWgtFun'));
h = p.Results.derivstep;
if isempty(inputtype) % either linear or nonlinear
    if isempty(betaCoef0) || isModelLinear(modelfun,betaCoef0,x,h)
        if isRobust
            lsrtype = 3;
        else 
            lsrtype = 1;
        end
    else
        if isRobust
            lsrtype = 4;
        else 
            lsrtype = 2;
        end
    end
elseif any(strcmp(inputtype,{'ols','wls','gls','lin','linear'}))
    lsrtype = 1;
elseif any(strcmp(inputtype,{'nlin','nonlinear'}))
    lsrtype = 2;    
elseif any(strcmp(inputtype,{'total','tls'}))
    lsrtype = 5;
elseif any(strcmp(inputtype,{'robust'})) %nonlinear and linear robust 
    if isempty(betaCoef0) || isModelLinear(modelfun,betaCoef0,x,h)
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

%% IsValid Functions
function isLinear = isModelLinear(modelfun,betacoef,x,h)
%% Compute the Numerical Hessian to determine linearity of modelfun
%check for each point
for i=1:size(x,1)
    Jfun = @(bn)(modelfun(bn,x(i,:)));
    J = @(b) (calcPartials(Jfun,b,h));
    H = calcPartials(J,betacoef',h); %calculate hessian
    isLinear = ~any(H(:));
    if ~isLinear
       break; 
    end
end
end

function [robustTune, robustWgtFun, isvalid] = getRobust(weightfun, tune)
%% the weightfunction can be either a string of a function
% determine which and make sure its valid
% return the tune and weightfun

robustTune = tune;
robustWgtFun = weightfun;

inputWgtFun = weightfun;
if ischar(inputWgtFun)
    isvalid = true;
    switch inputWgtFun
        case 'andrews'
            robustWgtFun = @(x) (abs(x)<pi) .* sin(x)./x;
            defaultTune = 1.339;
        case 'bisquare'
            robustWgtFun = @(x) (abs(x)<1) .* (1-x.^2).^2;
            defaultTune = 4.685;
        case 'cauchy'
            robustWgtFun = @(x) 1./(1+x.^2);
            defaultTune = 2.385;
        case 'fair'
            robustWgtFun = @(x) 1./(1+abs(x));
            defaultTune = 1.400;
        case 'huber'
            robustWgtFun = @(x) 1./(max(1,abs(x)));
            defaultTune = 1.345;
        case 'logistic'
            robustWgtFun = @(x) tanh(x)./x;
            defaultTune = 1.205;
        case 'talwar'
            robustWgtFun = @(x) double(abs(x)<1);
            defaultTune = 2.795;
        case 'welsch'
            robustWgtFun = @(x) exp(-1*x.^2);
            defaultTune = 2.985;
        otherwise
            isvalid = false; %don't know this function
    end
    if isempty(robustTune)
        robustTune = defaultTune;
    end
elseif ~isempty(weightfun)
 %user explicitly input a robust weight function
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

function [Sx, isinputcovariance] = getCovariance(lstype,weights,x)
% return a covariance depending on the type of input
% vector is "weights", where covariance is inverse of those on the diagonal
% matrix is assumed to be a covariance
nEqn = size(x,1);
nVars = numel(x);
if isempty(weights)
    if lstype == 5 % total least squares
        Sx = speye(nVars);
    else
        Sx = speye(nEqn);
    end
    isinputcovariance = false;
elseif isvector(weights) %matrix is weights
    Sx = diag(1./weights);
    isinputcovariance = false;
else %matrix is covariance
    Sx = weights;
    isinputcovariance = true;
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

nEqn = size(modelfun(b,x),1)/size(x,1);

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



function K = calcK(modelfun,betacoef,x,y)
% calculate the K matrix for each iteration
K = y - modelfun(betacoef,x);

end