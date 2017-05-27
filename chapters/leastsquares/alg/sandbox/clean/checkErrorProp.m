function checkErrorProp
%% SLOPOV Volume Estimation
modelfunVolume = @(x)(x(1)*x(2)*x(3));
x_numeric = [10.1 4.7 6.3];
Sxx_SLOPOV_numeric = [0.25 0.03 0.1].^2;

[covY,Y] = calcErrorProp(modelfunVolume,x_numeric,Sxx_SLOPOV_numeric);
%check
checkVals(Y,299.0610,1e-10);
checkVals(covY,80.97491446,1e-6);

%% SLOPOV Volume Estimation with Explicit Partial Function
modelfunVolume = @(x)(x(1)*x(2)*x(3));
partialfunYXVolume = @(x)([x(2)*x(3) x(1)*x(3) x(1)*x(2)]);
[covY,Y] = calcErrorProp(modelfunVolume,x_numeric,Sxx_SLOPOV_numeric,...
    'JacobianYX',partialfunYXVolume);
%check
checkVals(Y,299.0610,1e-10);
checkVals(covY,80.97491446,1e-10);

%% GLOPOV y=mx+b with modelfun(x)
modelfunMxplusb = @(x)([x(1)*x(3)+x(2);x(1)*x(4)+x(2)]);

Sxx_2 = [0.2 -1 0 0;-1 10 0 0;0 0 0 0;0 0 0 0.2];
x2 = [1.25 0.3 3 5];
[covY,Y] = calcErrorProp(modelfunMxplusb,x2,Sxx_2);
%check
checkVals(Y,[4.05 6.55]',1e-10);
checkVals(covY,[5.8 5;5 5.3125],1e-8);

%% GLOPOV y=mx+b with modelfun(x) and Explicit Partial X
modelfunMxplusb = @(x)([x(1)*x(3)+x(2);x(1)*x(4)+x(2)]);
partialfunYXmxplusb = @(x) [x(3) 1 x(1) 0;x(4) 1 0 x(1)];

Sxx_2 = [0.2 -1 0 0;-1 10 0 0;0 0 0 0;0 0 0 0.2];
x2 = [1.25 0.3 3 5];
[covY,Y] = calcErrorProp(modelfunMxplusb,x2,Sxx_2,...
    'JacobianYX',partialfunYXmxplusb);
%check
checkVals(Y,[4.05 6.55]',1e-10);
checkVals(covY,[5.8 5;5 5.3125],1e-8);

%% GLOPOV y=mx+b with modelfun(b,x)
modelfunLinear = @(b,x) b(1)*x + b(2);
b = [1.25 0.3]';
Sbb = [0.2 -1;-1 10];
x = [3 5]';
Sxx = [0 0;0 0.2];

[covY,Y] = calcErrorProp(modelfunLinear,x,Sxx,b,Sbb);
%check
checkVals(Y,[4.05 6.55]',1e-10);
checkVals(covY,[5.8 5;5 5.3125],1e-8);

%% GLOPOV y=mx+b with modelfun(b,x) and Explicit Partial B
modelfunLinear = @(b,x) b(1)*x + b(2);
partialfunJBLinear = @(b,x) [x ones(size(x))];

b = [1.25 0.3]';
Sbb = [0.2 -1;-1 10];
x = [3 5]';
Sxx = [0 0;0 0.2];

[covY,Y] = calcErrorProp(modelfunLinear,x,Sxx,b,Sbb,...
    'JacobianYB',partialfunJBLinear);
%check
checkVals(Y,[4.05 6.55]',1e-10);
checkVals(covY,[5.8 5;5 5.3125],1e-8);

%% GLOPOV y=mx+b with modelfun(b,x) and Explicit Partial X and Partial B
modelfunLinear = @(b,x) b(1)*x + b(2);
partialfunJBLinear = @(b,x) [x ones(size(x))];
partialfunJXLinear = @(b,x) diag(ones(size(x))*b(1));

b = [1.25 0.3]';
Sbb = [0.2 -1;-1 10];
x = [3 5]';
Sxx = [0 0;0 0.2];

[covY,Y] = calcErrorProp(modelfunLinear,x,Sxx,b,Sbb,...
    'JacobianYB',partialfunJBLinear,'JacobianYX',partialfunJXLinear);
%check
checkVals(Y,[4.05 6.55]',1e-10);
checkVals(covY,[5.8 5;5 5.3125],1e-10);

%%
fprintf('All Tests Passed\n');
end

function checkVals(x,y,thresh)
dx = x-y;
isbad = any(abs(dx(:))>thresh);
if isbad
   error('Check Failed'); 
end
end