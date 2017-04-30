function [betacoef,R,J,CovB,MSE,ErrorModelInfo] = lsr(x,y,modelfun,betacoef0,options)
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
%   - 'So2is1'               : false (default), true 
%   - 'beta0Covariance'      : n/a (default), (m x m) covariance

%% temporary hardcode logic
W = options;
maxiter = 100;
%% Calculate NonLinear Least Squares
betacoef = betacoef0;                         % set first guess at unknowns

So2 = inf; dSo2 = 1; iter = 0;  % initialize while loop

m = size(x,1);            % number of observations
n = numel(betacoef0);                  % number of unknowns
dof = m-n;                      % degrees of freedom
fprintf('iter :        So2                 X(1)                ... \n');
while dSo2>0 && iter<maxiter %loop until So2 increases or exceed 100 iteration
    J = calcJ(modelfun,betacoef,x);
    K = calcK(modelfun,betacoef,x,y);
    dbetacoef = (J'*W*J)\J'*W*K;       % Loop Delta Estimate
    betacoef = betacoef + dbetacoef;                 % Loop Estimate
    V = K;                      % Residuals
    dSo2 = So2 - V'*W*V/dof;    % Change in Reference Variance
    So2 = (V'*W*V)/dof;         % Reference Variance
    iter = iter + 1;

    % print status to screen
    fprintf('%3.0f : ',iter);
    fprintf('%20.10f', So2);
    fprintf('%20.10f', betacoef);
    fprintf('\n');
end
ErrorModelInfo.X = betacoef;                          % Unknowns
ErrorModelInfo.V = V;                          % Residuals
ErrorModelInfo.So2 = So2;                      % Reference Variance
ErrorModelInfo.iter = iter;                    % total iterations

ErrorModelInfo.Q = inv(J'*W*J);                % cofactor
ErrorModelInfo.Sx = So2 * ErrorModelInfo.Q;           % covariance of unknowns
ErrorModelInfo.Sl = J * ErrorModelInfo.Sx * J';       % covariance of observations
ErrorModelInfo.stdX = sqrt(diag(ErrorModelInfo.Sx));  % std of solved unknowns
ErrorModelInfo.Lhat = J * betacoef;                   % predicted L values
ErrorModelInfo.RMSE = sqrt(V'*V/m);            % RMSE

%% Handle Output Variables
betacoef = ErrorModelInfo.X;   % unknowns
CovB = ErrorModelInfo.Sx; % covariance of unknowns
R = V;
MSE = So2;
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

