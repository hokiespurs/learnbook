function [X,Sx,lsrinfo] = lsrlin(A,L,s)
% LSRLIN solve linear least squares regression with optional weights
%   Solve an overconstrained, linear least squares regression by inputting
%   the A(Beta) and L(Y) matrices and an optional stochastic model.  The
%   stochastic model can either be a vector consisting of a weight for each
%   observation, or a covariance matrix for the observations.
% 
% Inputs:
%   - A : mxn        : A Matrix
%   - L : mx1        : L Matrix
%   - s : mx1 or mxm : stochastic model
% 
% Outputs:
%   - X       : nx1 :  Calculated Unkowns
%   - Sx      : nxn :  Covariance of Unknowns
%   - lsrinfo :  1  :  structure with extra least squares regression info
% 
% Examples:
%   x = [0 1 2 3 4];
%   y = [5 1 7 13 24];
%   weights = [1 10 100 5 1];
%   A = [x(:).^2 x(:) ones(size(x(:)))];
%   L = y(:);
%   [X,Sx,lsrinfo] = lsrlin(A,L,weights); %perform weighted least squares
%
% Dependencies:
%   - n/a
% 
% Toolboxes Required:
%   - n/a
% 
% Author        : Richie Slocum
% Email         : richie@cormorantanalytics.com
% Date Created  : 29-Mar-2017
% Date Modified : 29-Mar-2017
%% Store inputs
lsrinfo.inputs.A = A;
lsrinfo.inputs.L = L;
if nargin==2
    lsrinfo.inputs.s = [];
else
    lsrinfo.inputs.s = s; %stochastic model
end
%% Determine how to handle optional W matrix
if nargin==2 %No stochastic input
    W = eye(numel(L)); %equal weight
else
   [wm,wn]=size(s);
   if wm==1 || wn==1 %Weight matrix input is a vector
       W = diag(s); %each observation is weighted based on vector
   else %Weight matrix input is covariance
       W = inv(s); %weight is inverse of covariance matrix
   end
end

%% Calculate Linear Least Squares
m = numel(L);                   % number of observations
n = size(A,2);                  % number of unknowns
dof = m-n;                      % degrees of freedom
X = (A'*W*A)\A'*W*L;            % unknowns
V = A * X - L;                  % residuals
So2 = V'*W*V/dof;               % Reference Variance
Q = inv(A'*W*A);                % cofactor
Sx = So2 * Q;                   % covariance of unknowns
Sl = A * Sx * A';               % covariance of observations
stdX = sqrt(diag(Sx));          % std of solved unknowns
Lhat = A * X;                   % predicted L values
r2 = var(Lhat)/var(L);          % R^2 Skill
RMSE = sqrt(V'*V/m);            % RMSE

%% Populate lsrinfo Variables (didnt do above for clarity)
lsrinfo.m=m;
lsrinfo.n=n;
lsrinfo.dof=dof;
lsrinfo.X=X;
lsrinfo.V=V;
lsrinfo.So2=So2;
lsrinfo.Q=Q;
lsrinfo.Sx=Sx;
lsrinfo.Sl=Sl;
lsrinfo.stdX=stdX;
lsrinfo.Lhat=Lhat;
lsrinfo.r2=r2;
lsrinfo.RMSE=RMSE;

end