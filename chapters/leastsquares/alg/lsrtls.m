function [X,Sx,lsrinfo] = lsrtls(Jfun,Kfun,Bfun,Xo,s)
%% LSRNLIN solves linear least squares adjustment with optional weights
%
% Input
% Jfun : jacobian
% Kfun : l-F(x)
% Bfun : Jacobian wrt residuals
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
lsrinfo.inputs.Bfun = Bfun;
lsrinfo.inputs.Xo = Xo;
if nargin==4
    lsrinfo.inputs.s = [];
else
    lsrinfo.inputs.s = s;
end

%% Determine how to handle optional W matrix (Make it covariance)
if nargin==4 %No stochastic input
    S = eye(size(Bfun(Xo),2)); %equal weight case
else
   [wm,wn]=size(s);
   if wm==1 || wn==1 %Weight matrix input is a vector
       S = inv(s); %inverse weight (another inverse happens later)
   else %Weight matrix input is covariance
       S = s; %weight is covariance matrix (inverse happens later)
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
    B = Bfun(X);
    J = Jfun(X);
    K = Kfun(X);
    Weq = inv(B*S*B');                % equivalent weight matrix
    dX = (J'*Weq*J)\J'*Weq*K;         % Loop Delta Estimate
    X = X + dX;                       % Loop Estimate
    Veq = K;                          % Residuals
    dSo2 = So2 - Veq'*Weq*Veq/dof;    % Change in Reference Variance
    So2 = (Veq'*Weq*Veq)/dof;         % Reference Variance
    iter = iter + 1;                  % Count Iterations
    % print status to screen
    fprintf('%3.0f : ',iter);
    fprintf('%20.10f',So2);
    fprintf('%20.10f',X);
    fprintf('\n');
end
lsrinfo.X = X;                          % Unknowns
lsrinfo.Veq = Veq;                      % Equivalent Residuals
lsrinfo.Vobs = S * B' * Weq * Veq;      % Observation Residuals
lsrinfo.So2 = So2;                      % Reference Variance
lsrinfo.iter = iter;                    % total iterations

lsrinfo.Q = inv(J'*Weq*J);              % cofactor
lsrinfo.Sx = So2 * lsrinfo.Q;           % covariance of unknowns
lsrinfo.Sl = J * lsrinfo.Sx * J';       % covariance of observations
lsrinfo.stdX = sqrt(diag(lsrinfo.Sx));  % std of solved unknowns
lsrinfo.Lhat = J * X;                   % predicted L values
lsrinfo.RMSE = sqrt(Veq'*Veq/m);        % RMSE

%% Handle Output Variables
X = lsrinfo.X;   % unknowns
Sx = lsrinfo.Sx; % covariance of unknowns
end