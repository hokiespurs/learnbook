function [Syy,y] = calcErrorProp(modelfun,x,Sxx,varargin)
% input options
%
% calcErrorProp(modelfun,x,Sxx)
% calcErrorProp(modelfun,x,Sxx,b,Sbb)
%
% Options
% - derivStep
% - partialfunJyb
% - partialfunJyx
%  
% modelfun is the function handle which takes x as an input @(x)
% x is data
% Sxx is covariance of x data
% b is beta coef
% Sbb is covariance of b data
% partialFun is analytical partial derivative

%% Input Parsing/Checks
narginchk(3,11);
nobs = numel(x);

% default values
defaultb   = [];
defaultSbb = [];
defaultDerivStep  = eps^(1/3);
defaultJacobianYB = [];
defaultJacobianYX = [];

% Input check functions
checkModelfun      = @(X) isa(X,'function_handle');
checkX             = @(X) isnumeric(X) && all(~isnan(X(:)));
checkSxx           = @(X) isnumeric(X) && all(~isnan(X(:))) ...
    && (isCovariance(X,nobs) || (isvector(X) && numel(X)==nobs));
checkB             = @(X) isnumeric(X) && all(~isnan(X(:)));
checkSbb           = @(X) isnumeric(X) && all(~isnan(X(:))) ...
    && (isvector(X) || isCovariance(X)); %no info about size yet
checkDerivStep     = @(X) isnumeric(X) && isscalar(X) && X>0;
checkJacobianYB    = @(X) isa(X,'function_handle');
checkJacobianYX    = @(X) isa(X,'function_handle'); 

% Parse Inputs
p = inputParser;

addRequired(p,'modelfun' ,checkModelfun);
addRequired(p,'x'        ,checkX);
addRequired(p,'Sxx'      ,checkSxx);

addOptional(p,'b'   ,defaultb   ,checkB);
addOptional(p,'Sbb' ,defaultSbb ,checkSbb);

addParameter(p,'DerivStep'  ,defaultDerivStep  ,checkDerivStep);
addParameter(p,'JacobianYB' ,defaultJacobianYB ,checkJacobianYB);
addParameter(p,'JacobianYX' ,defaultJacobianYX ,checkJacobianYX);

parse(p,modelfun,x,Sxx,varargin{:});

%% Get variables out of structure
b = p.Results.b(:);
Sbb = p.Results.Sbb;
h = p.Results.DerivStep;

% if no covariance for b then make it zeros
if ~isempty(b) && isempty(Sbb)
    Sbb = sparse(zeros(numel(p.Results.b)));    
end

% turn vector Sbb and Sxx into covariances
if isvector(Sbb)
   Sbb = spdiag(Sbb); 
end
if isvector(Sxx)
    Sxx = spdiag(Sxx);
end

if isempty(Sbb) % no beta into model
    if isempty(p.Results.JacobianYX)
        Jyx = calcPartials(modelfun,x,h);
        Jyb = [];
    else
        Jyx = p.Results.JacobianYX(x);
        Jyb = [];
    end
else
    if isempty(p.Results.JacobianYX)
        Jyx = calcJYX(modelfun,b,x,h);
    else
        Jyx = p.Results.JacobianYX(b,x);
    end
    if isempty(p.Results.JacobianYB)
        Jyb = calcJYB(modelfun,b,x,h);
    else
        Jyb = p.Results.JacobianYB(b,x);
    end    
end

%% Check Getting Y data
try
   if isempty(Sbb)
       y = modelfun(x);
   else
       y = modelfun(b,x);
   end
catch
    error('Error calling modelfun');
end

%% Do GLOPOV
Jybyx = [sparse(Jyb) sparse(Jyx)];
Covbx = blkdiag(sparse(Sbb),sparse(Sxx));
Syy = Jybyx * Covbx * Jybyx';

end

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

function Y = spdiag(X)
% make a sparse matrix on the diagonal using vector X
n = numel(X);
Y = sparse(1:n,1:n,X);
end

function isvalid = isCovariance(C,ndim)
% check if a covariance matrix is valid
if nargin==1
   ndim = size(C,1); 
end
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