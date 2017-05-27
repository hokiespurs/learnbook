function exampleErrorProp
%% SLOPOV Volume Estimation
modelfunVolume = @(x)(x(1)*x(2)*x(3));
x_numeric = [10.1 4.7 6.3];
Sxx_SLOPOV_numeric = [0.25 0.03 0.1].^2;

[covY,Y] = calcErrorProp(modelfunVolume,x_numeric,Sxx_SLOPOV_numeric);

%% GLOPOV y=mx+b with modelfun(b,x)
modelfunLinear = @(b,x) b(1)*x + b(2);
b = [1.25 0.3]';
Sbb = [0.2 -1;-1 10];
x = [3 5]';
Sxx = [0 0;0 0.2];

[covY,Y] = calcErrorProp(modelfunLinear,x,Sxx,b,Sbb);

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

%% Calculate Position, Velocity and Acceleration of Projectile
% initial velocity, time, acceleration
projectilefun = @(x)[x(1)*x(2)+0.5*x(3)*x(2)^2;x(1)+x(3)*x(2)];
x = [1 2 3];
sx = [0.1 0.3 2];
[covY,Y] = calcErrorProp(projectilefun,x,sx);

%% Calculating Photogrammetric Projection of Pixels onto Plane

end
