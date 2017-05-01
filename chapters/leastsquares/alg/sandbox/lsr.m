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
%% PARSE INPUTS
p = inputParser;
% default values
defaultType = 'nonlinear';
defaultStochastic = [];
defaultJacobian = @(b,x) calcJ(modelfun,b,x);
defaultBfunction = @(b,x) calcB(modelfun,b,x);
defaultnoscale = false;
defaultBeta0Cov = [];
defaultVerbose = false;
defaultMaxIter = 100;
defaultBetacoefAsObs = false;
defaultchi2 = 0.05;
% expected options
expectedType = {'ols','wls','gls','lin','linear',...
    'nlin','nonlinear',...
    'total','tls',...
    'robust'};
% add input validation
isfunction = @(x)(isa(x,'function_handle'));

addRequired(p,'x',@isnumeric);
addRequired(p,'y',@(x)(isnumeric(x) & size(x,2)==1));
addRequired(p,'modelfun',isfunction);

addOptional(p,'betacoef0',defaultBeta0Cov,@isnumeric);

addParameter(p,'type',defaultType,@(x) any(validatestring(x,expectedType)));
addParameter(p,'stochastic',defaultStochastic,@isnumeric)
addParameter(p,'analyticalJ',defaultJacobian,isfunction)
addParameter(p,'analyticalB',defaultBfunction,isfunction)
addParameter(p,'noscale',defaultnoscale,@islogical)
addParameter(p,'beta0covariance',defaultBeta0Cov,@isnumeric)
addParameter(p,'verbose',defaultVerbose);
addParameter(p,'maxiter',defaultMaxIter,@(x)(isscalar(x) & isnumeric(x)));
addParameter(p,'betacoef0asObs',defaultBetacoefAsObs,@islogical);
addParameter(p,'chi2alpha',defaultchi2,@(x)(isscalar(x) & isnumeric(x)));
parse(p,x,y,modelfun,varargin{:});
ErrorModelInfo.parser = p; % include function input in metadata output

%% ADDITIONAL INPUT CHECKS
% Determine Level of output to screen
isverbose = p.Results.verbose;

% check what type of least squares
[ErrorModelInfo.meta.type, islinear, isnonlinear, istls] = getLsrType(p,expectedType,isverbose);

% determine nObsEqns, nObservations, nPredictors, nEqnPerObservations
nObsEqns = numel(y);
nPredictors = numel(x);
nBetacoef = calcNbetacoef(modelfun,x);
dof = nObsEqns-nBetacoef; %degrees of freedom
if dof<=0
   error('Must be an overconstrained system of equations'); 
end
% assemble stochastic model into covariance matrix covY (or covX for TLS)
[ErrorModelInfo.meta.stochastic,covX,covY] = calcStochastic(p,isverbose,nPredictors,nObsEqns,istls);

% do checks specific to nonlinear and tls functions
if (isnonlinear || istls)
    if isempty(p.Results.betacoef0) %ensure initial beta guess is  for nonlinear and tls
        error('Need an initial guess at Beta coefficients for Nonlinear or TLS');
    end
   betacoef0 = p.Results.betacoef0(:); %ensure one column
   [ErrorModelInfo.meta.Jtype,Jfun]=checkJ(p,betacoef0,x,y,isverbose);
   if istls
       [ErrorModelInfo.meta.Btype,Bfun]=checkB(p,isverbose,betacoef0,x,y);
   end
end

% check and populate Beta Covariance
betacoef0asobs = p.Results.betacoef0asObs;
betacoefcov = p.Results.beta0covariance;
[betacoef0asobs,betacoefcov] = checkBetacoefCov(p,betacoef0asobs,betacoefcov,nBetacoef,islinear);

%% DO LEAST SQUARES
if islinear
%% Calculate Linear Least Squares
A = calcJ(modelfun,zeros(nBetacoef,1),x);
L = y;

if betacoef0asobs %add betacoef0 as observation equations
    if isverbose
       fprintf('Using betacoef0 as observation equations: yes\n'); 
    end
    
    L(end+1:end+nBetacoef)=p.Results.betacoef0;
    A = [A;eye(nBetacoef)];
    covY = blkdiag(covY,betacoefcov);
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
    fprintf('So2=%.4f\n',So2);
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
ErrorModelInfo.Sl = A * ErrorModelInfo.covB * A'; % covariance of observations
ErrorModelInfo.stdX = sqrt(diag(ErrorModelInfo.covB));  % std of solved unknowns
ErrorModelInfo.Lhat = A * betacoef;                   % predicted L values
ErrorModelInfo.RMSE = sqrt(V'*V/nObsEqns);        % RMSE
ErrorModelInfo.r2 = var(ErrorModelInfo.Lhat)/var(L);  % R^2 Skill
% handle output variables
R = V;
MSE = So2;
CovB = ErrorModelInfo.covB;
J = A;
elseif isnonlinear || istls
%% Calculate NonLinear Least Squares and TLS
    maxiter = p.Results.maxiter;
    betacoef = betacoef0;                  % set first guess at unknowns
    if isverbose
        if betacoef0asobs
       fprintf('Using betacoef0 as observation equations: yes\n'); 
        else
       fprintf('Using betacoef0 as observation equations: no\n'); 
        end        
    end
    %initialize while loop parameters
    So2 = inf;
    dSo2 = 1;
    iter = 0;
    if isverbose
        fprintf('\niter :        So2        '); %
        fprintf('         betacoef(%.0f)',1:nBetacoef);
        fprintf('\n');
    end
    while dSo2>0 && iter<maxiter %loop until So2 increases or exceed 100 iteration
            J = Jfun(betacoef,x);
            K = calcK(modelfun,betacoef,x,y);
            if istls
                B = Bfun(betacoef,x);
                W = inv(B*covX*B');            %equivalent weight matrix
            else
                W = inv(covY);
            end
            
            if betacoef0asobs
                if istls
                    [J,K,W,B]=addBetacoef0TLS(J,K,B,covX,betacoef0,betacoefcov,betacoef);
                else
                    [J,K,W]=addBetacoef0Nlin(J,K,covY,betacoef0,betacoefcov,betacoef);
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
    ErrorModelInfo.betacoef = betacoef;
    ErrorModelInfo.R = V;
    if istls
        if betacoef0asobs
            covX = blkdiag(covX,betacoefcov);
        end
       ErrorModelInfo.Robs = covX * B' * W * V; 
    end
    ErrorModelInfo.MSE = So2;
    ErrorModelInfo.iter = iter;
    if isverbose
        fprintf('So2=%.4f\n',So2);
    end
    ErrorModelInfo.Q = inv(J'*W*J);              % cofactor
    if p.Results.noscale %force it not to scale covariance
       if isverbose
           fprintf('Covariance is not scaled\n');
       end
       ErrorModelInfo.covB = inv(J'*W*J);
    else
       if isverbose
           fprintf('Covariance is scaled\n');
       end
       ErrorModelInfo.covB = So2*inv(J'*W*J);
    end
    ErrorModelInfo.Sl = J * ErrorModelInfo.covB * J'; % covariance of observations
    ErrorModelInfo.stdX = sqrt(diag(ErrorModelInfo.covB));  % std of solved unknowns
    ErrorModelInfo.Lhat = J * betacoef;                   % predicted L values
    ErrorModelInfo.RMSE = sqrt(V'*V/nObsEqns);        % RMSE
    % handle output variables
    R = V;
    MSE = So2;
    CovB = ErrorModelInfo.covB;
else %robust
%% Robust Least Squares    

end
%% Do Chi2 Test
if strcmp(ErrorModelInfo.meta.stochastic, 'Covariance')
    chi2 = dof * So2;
    alphaval = p.Results.chi2alpha;
    chi2low = chi2inv(alphaval/2,dof);
    chi2high = chi2inv(1-alphaval/2,dof);
    if chi2>chi2low && chi2<chi2high
        if p.Results.noscale == false
           if isverbose
               warning('consider setting "noscale" equal to true, because the chi2 test passed');
           end
        end
        if isverbose
            fprintf('Chi2 Test Passed, We cant say So2 is not equal 1 ');
        end
        ErrorModelInfo.chi2.testpass = true;
        ErrorModelInfo.chi2.computedchi2 = chi2;
        ErrorModelInfo.chi2.chi2low = chi2low;
        ErrorModelInfo.chi2.chi2high = chi2high;
    else
        if isverbose
            if chi2>chi2high
                warning('Based on the chi2 test, the stochastic model is either incorrect (overly optimistic), or the data likely contains outliers\n');
            else
                warning('Based on the chi2 test, the stochastic model is overestimating the errors.  Your measurement accuracy is better than the stochastic mode');
            end
            fprintf('Chi2 Test Failed, So2 is not equal to 1 ');
        end
        ErrorModelInfo.chi2.testpass = false;
        ErrorModelInfo.chi2.computedchi2 = chi2;
        ErrorModelInfo.chi2.chi2low = chi2low;
        ErrorModelInfo.chi2.chi2high = chi2high;
    end
    if isverbose
        fprintf('at %.2f level of confidence\n',alphaval);
    end
else
    if isverbose
        fprintf('Stochastic model isnt covariances, so chi2 test isnt valid\n');
    end    
end

end

function nBetacoef = calcNbetacoef(modelfun,x)
   for iTestBeta=1:100 %loop and test different nBetacoef and see what works
       try
           modelfun(zeros(iTestBeta,1),x(2,:));
           nBetacoef = iTestBeta;
           break;
       catch
           % keep looping until modelfun doesnt throw an error
           % maybe I shouldnt have tried to force linear into this function
       end
   end
   if nBetacoef ==100
      error('Unable to determine the number of beta coeficients'); 
   end
end

function [J,K,W,B]=addBetacoef0TLS(J,K,B,covX,betacoef0,betacoefcov,betacoef)

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

function J = calcJ(modelfun,betacoef,x)
% calculate the Jacobian for each iteration 
% Jacobian is the partials wrt the betacoef regression coefficients
h = eps^(1/3); % optimal for central difference
Jfun = @(bn)(modelfun(bn,x));
J = calcPartials(Jfun,betacoef',h);

end

function B = calcB(modelfun,betacoef,x)
% calculate the B matrix for each iteration (only TLS)
% B matrix is partials wrt predictor variables
h = eps^(1/3); % optimal for central difference
Bfun = @(xn)(modelfun(betacoef,xn));
B = calcPartials(Bfun,x,h);

nEqn = size(modelfun(betacoef,x),1)/size(x,1);
B = bumphdiag(B,nEqn);
end

function K = calcK(modelfun,betacoef,x,y)
% calculate the K matrix for each iteration
K = y - modelfun(betacoef,x);

end

function [stochastictype,covX,covY] = calcStochastic(p,isverbose,nPredictors,nObsEqns,istls)
%one of these won't get called, and doesnt get used later on
covX = nan;
covY = nan;

if isempty(p.Results.stochastic) %default to identity matrix
    stochastictype = 'N/A';
    if isverbose
        fprintf('Stochastic Model : N/A\n');
    end
    if istls
        covX = eye(nPredictors);
    else
        covY = eye(nObsEqns);
    end
elseif isvector(p.Results.stochastic) %weight for each observation equation
    stochastictype = 'Weight Vector';
    if isverbose
       fprintf('Stochastic Model : Weight Vector\n'); 
    end
    if istls
        error('stochastic input must be a covariance matrix when doing tls, cant be a vector');
    end
    if numel(p.Results.stochastic)~=nObsEqns
        error('stochastic (weight) vector must be the same number of elements as y')
    end
    covY = inv(diag(p.Results.stochastic));
else % user input is covariance
    stochastictype = 'Covariance';
    if isverbose
       fprintf('Stochastic Model : Covariance\n'); 
    end
    if istls
        covX = p.Results.stochastic;
        if size(covX,1)~=nPredictors || ~issymmetric(covX)
            error('Invalid Covariance matrix, must be symmetric with size [%.0f x %.0f]',...
                nPredictors,nPredictors);
        end
    else
        covY = p.Results.stochastic;
        if size(covY,1)~=nObsEqns || ~issymmetric(covY)
            error('Invalid Covariance matrix, must be symmetric with size [%.0f x %.0f]',...
                nPredictors,nPredictors);
        end
    end
end

end

function [Jtype,Jfun]=checkJ(p,betacoef0,x,y,isverbose)
    if ismember('analyticalJ',p.UsingDefaults)
        Jtype = 'numerical';
        if isverbose
            fprintf('Jacobian Function: Numerical\n');
        end
    else
        Jtype = 'analytical';
        if isverbose
            fprintf('Jacobian Function: Analytical\n');
        end
    end
    Jfun = p.Results.analyticalJ;
    % make sure the function works
    try
       foo=Jfun(betacoef0,x)*betacoef0;
       if any(size(foo)~=size(y))
           error('J*beta should return the same size as the response(y)');
       end
    catch Jerr
        if ismember('analyticalJ',p.UsingDefaults) %default value
            fprintf('Calculating Jacobian of Modelfun throws an error\n');
        else
            fprintf('J Function throws an error\n');
        end
        rethrow(Jerr);
    end
end

function [lsrtype, islinear, isnonlinear, istls]=getLsrType(p,expectedType,isverbose)
hline = repmat('-',1,50);
if ismember(p.Results.type,expectedType(1:5))
    lsrtype = 'linear';
    if isverbose
       fprintf('\n%s\nPerforming Linear Least Squares\n%s\n',hline,hline); 
    end
    islinear = true;
    isnonlinear = false;
    istls = false;
elseif ismember(p.Results.type,expectedType(6:7))
    lsrtype = 'nonlinear';
    if isverbose
       fprintf('\n%s\nPerforming NonLinear Least Squares\n%s\n',hline,hline); 
    end
    islinear = false;
    isnonlinear = true;
    istls = false;
else
    lsrtype = 'total';
    if isverbose
       fprintf('\n%s\nPerforming Total Least Squares\n%s\n',hline,hline); 
    end
    islinear = false;
    isnonlinear = false;
    istls = true;
end
end

function [Btype,Bfun]=checkB(p,isverbose,betacoef0,x,y)
if ismember('analyticalB',p.UsingDefaults)
    Btype = 'numerical';
    if isverbose
        fprintf('B Function       : Numerical\n');
    end
else
    Btype = 'analytical';
    if isverbose
        fprintf('B Function       : Analytical\n');
    end
end
Bfun = p.Results.analyticalB;
% make sure the function works
try
    foo=Bfun(betacoef0,x)*x(:);
    if sum(size(foo)~=size(y))>0
        error('B function*x(:) should return the same size as the response(y)');
    end
catch Berr
    if ismember('analyticalB',p.UsingDefaults) %default value
        fprintf('Calculating B function of Modelfun throws an error\n');
    else
        fprintf('B Function throws an error\n');
    end
    rethrow(Berr);
end
end

function [betacoef0asobs,betacoefcov] = checkBetacoefCov(p,betacoef0asobs,betacoefcov,nBetacoef,islinear)
if ~isempty(betacoefcov)
    betacoef0asobs = true;
    if isempty(p.Results.betacoef0) % covariance but no initial guess
        error('Need to input beta coefficients to add as observation equations with the input covariance');
    end
    if ~issymmetric(betacoefcov) || size(betacoefcov,1)~=nBetacoef
        error('beta0covariance must be a valid covariance matrix with size [%.0f x %.0f)',...
            nBetacoef,nBetacoef);
    end
else
    betacoefcov = eye(nBetacoef);
end
% Checks for betacoef0 as observations
if betacoef0asobs && isempty(betacoefcov) && ~isempty(p.Results.stochastic)
   error('When stochastic input, betacoefasObs = true, need betacoef0covariance');
elseif betacoef0asobs && ~isempty(betacoefcov) && isempty(p.Results.stochastic)
   error('When stochastic not input, betacoefasObs = true, cant have betacoef0covariance because itd be weighting everything else as 1, and the beta0covariacne differently');
elseif islinear && isempty(p.Results.betacoef0) && betacoef0asobs
    error('Need to input betacoef0 in order to use it as an observation');
elseif islinear && ~betacoef0asobs && ~isempty(p.Results.betacoef0)
    warning('betacoef0 isnt doing anything for this linear least squares solution, use betacoef0asobs to include it');
end
end