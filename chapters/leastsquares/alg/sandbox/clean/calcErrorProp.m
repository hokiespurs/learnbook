function [Y,covY] = calcErrorProp(modelfun,varargin)
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
% S is covariance of x data
% partialFun is analytical partial derivative

h = eps^(1/3);
if nargin==3
   Jyx = calcPartials(modelfun,x,h);
else
   Jyx = partialfun(x);
end
if isvector(Sxx)
    Sxx = diag(Sxx);
end
Y = modelfun(x);
covY = Jyx*Sxx*Jyx.';

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