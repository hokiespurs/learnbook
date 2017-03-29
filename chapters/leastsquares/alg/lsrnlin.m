function [X,Sx,lsrinfo] = lsrnlin(Jfun,Kfun,Xo,s)
%% LSRNLIN solves linear least squares adjustment with optional weights
%
% Input
% Jfun : jacobian
% Kfun : l-F(x)
% Xo : initial guess
% s: optional stochastic model: vector of weights or matrix of covariances
%
% Output
% X
% Sx
% lsainfo : extra output parameters
%% inputs
lsrinfo.inputs.Jfun = Jfun;
lsrinfo.inputs.Kfun = Kfun;
lsrinfo.inputs.Xo = Xo;
if nargin==3
    lsrinfo.inputs.s = [];
else
    lsrinfo.inputs.s = s;
end

%% Determine how to handle optional W matrix
if nargin==3 %No stochastic input
    W = eye(size(Kfun(Xo),1)); %equal weight
else
   [wm,wn]=size(s);
   if wm==1 || wn==1 %Weight matrix input is a vector
       W = diag(s); %each observation is weighted based on vector
   else %Weight matrix input is covariance
       W = inv(s); %weight is inverse of covariance matrix
   end
end

%% Calculate NonLinear Least Squares
X = Xo;                         % set first guess at unknowns

So2 = inf; dSo2 = 1; iter = 0;  % initialize while loop

m = size(Jfun(X),1);            % number of observations
n = numel(Xo);                  % number of unknowns
dof = m-n;                      % degrees of freedom
fprintf('iter :        So2                 X(1)                ... \n');

while dSo2>0 && iter<100 %loop until So2 increases or exceed 100 iteration
    J = Jfun(X);
    K = Kfun(X);
    dX = (J'*W*J)\J'*W*K;               % Loop Delta Estimate
    X = X + dX;                         % Loop Estimate
    V = K;                              % Loop Residuals
    dSo2 = So2 - (V'*W*V)/dof;          % Change in Reference Variance
    So2 = (V'*W*V)/dof;                 % Reference Variance
    iter = iter + 1;                    % Count Iterations
    % print status to screen
    fprintf('%3.0f : ',iter);
    fprintf('%20.10f',So2);
    fprintf('%20.10f',X);
    fprintf('\n');
end
lsrinfo.X = X;                          % Unknowns
lsrinfo.V = V;                          % Residuals
lsrinfo.So2 = So2;                      % Reference Variance
lsrinfo.iter = iter;                    % total iterations

lsrinfo.Q = inv(J'*W*J);                % cofactor
lsrinfo.Sx = So2 * lsrinfo.Q;           % covariance of unknowns
lsrinfo.Sl = J * lsrinfo.Sx * J';       % covariance of observations
lsrinfo.stdX = sqrt(diag(lsrinfo.Sx));  % std of solved unknowns
lsrinfo.Lhat = J * X;                   % predicted L values
lsrinfo.RMSE = sqrt(V'*V/m);            % RMSE

%% Handle Output Variables
X = lsrinfo.X;   % unknowns
Sx = lsrinfo.Sx; % covariance of unknowns
end