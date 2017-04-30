function [betacoef,R,J,CovB,MSE,ErrorModelInfo] = lsr(x,y,modelfun,betacoef0,varargin)
% LSR 
% Inputs:
%   - x              : Predictor variables
%   - y              : Response values
%   - modelfun       : Model function handle @modelfun(betacoef,X)
%   - betacoef0          : Initial coefficient values
%   - options        : Structure with optional parameters
% 
% Outputs:
%   - betacoef           : Estimated regression coefficients
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
% expected options
expectedType = {'ols','wls','gls','lin','linear',...
    'nlin','nonlinear',...
    'total','tls'};
% add input validation
isfunction = @(x)(isa(x,'function_handle'));

addRequired(p,'x',@isnumeric);
addRequired(p,'y',@isnumeric);
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

parse(p,x,y,modelfun,varargin{:});
%% Check Input Logic
%% Determine Level of output to screen
isverbose = p.Results.verbose;
hline = repmat('-',1,50);
%% check what type of least squares
if ismember(p.Results.type,expectedType(1:5))
    ErrorModelInfo.meta.type = 'linear';
    if isverbose
       fprintf('\n%s\nPerforming Linear Least Squares\n%s\n',hline,hline); 
    end
    islinear = true;
    isnonlinear = false;
    istls = false;
elseif ismember(p.Results.type,expectedType(6:7))
    ErrorModelInfo.meta.type = 'nonlinear';
    if isverbose
       fprintf('\n%s\nPerforming NonLinear Least Squares\n%s\n',hline,hline); 
    end
    islinear = false;
    isnonlinear = true;
    istls = false;
else
    ErrorModelInfo.meta.type = 'total';
    if isverbose
       fprintf('\n%s\nPerforming Total Least Squares\n%s\n',hline,hline); 
    end
    islinear = false;
    isnonlinear = false;
    istls = true;
end
%% Make sure y matrix is ok
if size(y,2)~=1
   error('Response vector y must be a (mx1) vector'); 
end
%% determine nObsEqns, nObservations, nPredictors, nEqnPerObservations
nObsEqns = numel(y);
nObservations = size(x,1);
nPredictors = numel(x);
nEqnPerObservations = nObsEqns/nObservations;
if islinear %no explicit info about betacoef is input for linear functions 
   for iTestBeta=1:100 %loop and test different nBetacoef and see what works
       try
           modelfun(zeros(iTestBeta,1),x);
           nBetacoef = iTestBeta;
           break;
       catch
           % keep looping until is works
           % maybe I shouldnt have tried to force linear into this function
       end
   end
else
   nBetacoef = numel(p.Results.betacoef0);
end
dof = nObsEqns-nBetacoef; %degrees of freedom
if dof<=0
   error('Must be an overconstrained system of equations'); 
end
%% assemble stochastic model into covariance matrix covY (or covX for TLS)
if isempty(p.Results.stochastic) %default to identity matrix
    ErrorModelInfo.meta.stochastic = 'N/A';
    if isverbose
        fprintf('Stochastic Model: N/A\n');
    end
    if istls
        covX = eye(nPredictors);
    else
        covY = eye(nObsEqns);
    end
elseif isvector(p.Results.stochastic) %weight for each observation equation
    ErrorModelInfo.meta.stochastic = 'Weight Vector';
    if isverbose
       fprintf('Stochastic Model: Weight Vector\n'); 
    end
    if istls
        error('stochastic input must be a covariance matrix when doing tls, cant be a vector');
    end
    if numel(p.Results.stochastic)~=nObsEqns
        error('stochastic (weight) vector must be the same number of elements as y')
    end
    covY = inv(diag(p.Results.stochastic));
else % user input is covariance
    ErrorModelInfo.meta.stochastic = 'Covariance';
    if isverbose
       fprintf('Stochastic Model: Covariance\n'); 
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
%% ensure initial beta guess is input
if (isnonlinear || istls) && isempty(p.Results.betacoef0)
    error('Need an initial guess at Beta coeffiicents for Nonlinear or TLS');
end
if isnonlinear || istls
   betacoef0 = p.Results.betacoef0(:); %ensure one column
end
%% test Jacobian
if isnonlinear || istls
    if ismember('analyticalJ',p.UsingDefaults)
        ErrorModelInfo.meta.Jtype = 'numerical';
        if isverbose
            fprintf('Jacobian Function: Numerical\n');
        end
    else
        ErrorModelInfo.meta.Jtype = 'analytical';
        if isverbose
            fprintf('Jacobian Function: Analytical\n');
        end
    end
    Jfun = p.Results.analyticalJ;
    % make sure the function works
    try
       foo=Jfun(betacoef0,x)*betacoef0;
       if sum(size(foo)~=size(y))>0
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
%% test B Matrix (only TLS)
if istls
    if ismember('analyticalB',p.UsingDefaults)
        ErrorModelInfo.meta.Btype = 'numerical';
        if isverbose
            fprintf('B Function: Numerical\n');
        end
    else
        ErrorModelInfo.meta.Btype = 'analytical';
        if isverbose
            fprintf('B Function: Analytical\n');
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
%% test Beta Covariance
betacoefcov = p.Results.beta0covariance;
if ~issymmetric(betacoefcov) || size(betacoefcov,1)~=numel(betacoef0)
    error('beta0covariance must be a valid covariance matrix with size [%.0f x %.0f)',...
        nBetacoef,nBetacoef);
end
%% DO LEAST SQUARES
if islinear
%% Calculate Linear Least Squares

else
%% Calculate NonLinear Least Squares and TLS
    maxiter = p.Results.maxiter;
    betacoef = betacoef0;                  % set first guess at unknowns

    %initialize while loop parameters
    So2 = inf;
    dSo2 = 1;
    iter = 0;
    if isverbose
        fprintf('iter :        So2        '); %
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
       ErrorModelInfo.Robs = covX * B' * W * V; 
    end
    ErrorModelInfo.MSE = So2;
    ErrorModelInfo.iter = iter;
    
    ErrorModelInfo.Q = inv(J'*W*J);              % cofactor
    if p.Results.noscale
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
end

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

