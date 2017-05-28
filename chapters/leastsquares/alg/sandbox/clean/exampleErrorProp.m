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
%exterior orientation [rpy(radians), camera xyz]
eo = [0.184250,3.107580,0.345273,-3.781813,-0.566171,100.193073];

%Exterior Orientation covariance
eocov = [0.000183,-0.000153,-0.000049,-0.021162,0.010974,-0.000909;...
-0.000153,0.000438,0.000073,0.051112,0.006242,-0.000299;...
-0.000049,0.000073,0.000023,0.008923,-0.001257,0.000111;...
-0.021162,0.051112,0.008923,6.063994,0.392444,0.000456;...
0.010974,0.006242,-0.001257,0.392444,1.436837,-0.097038;...
-0.000909,-0.000299,0.000111,0.000456,-0.097038,0.084431];

% Camera K matrix 
f  = 320; varf = 0;
cx = 320; varcx = 0;
cy = 240; varcy = 0;

% DSM Z value
dsmz = 0.25;vardsmz = 1;

% u pix
[upix,vpix] = meshgrid(1:10:640,1:50:480);

%uv covariance
uvcov = diag(0.5*ones(numel(upix)*2,1));

b = [eo f cx cy dsmz];
Sbb = blkdiag(eocov,varf,varcx,varcy,vardsmz);

x = [upix(:) vpix(:)];
Sxx = uvcov;

[covY,Y] = calcErrorProp(@OrthoToPlane,x,Sxx,b,Sbb);
Xw = Y(1:4:end);
Yw = Y(2:4:end);
Zw = Y(3:4:end);
Sw = Y(4:4:end);

sXYZS = diag(covY);
sX = sXYZS(1:4:end);
sY = sXYZS(2:4:end);
sZ = sXYZS(3:4:end);
sS = sXYZS(4:4:end);
sMagXY = sqrt(sX.^2+sY.^2);
figure(2);clf
scatter(Xw,Yw,10,sMagXY,'filled');axis equal
figure(3);
pcolor(covY);shading flat;set(gca,'ydir','reverse');caxis([-5 5]);
axis equal
end

function y = OrthoToPlane(b,x)
% parse out values into readable values
roll = b(1);
pitch = b(2);
yaw = b(3);
Xc = b(4);
Yc = b(5);
Zc = b(6);
f = b(7);
cx = b(8);
cy = b(9);
z0 = b(10);

u = x(:,1);
v = x(:,2);

% Solve for most probably Value
Xw=(Xc.*cx*sin(pitch) - Xc*u*sin(pitch) + f*z0*sin(roll)*sin(yaw) - v*z0*cos(roll)*sin(yaw) + Xc*f*cos(pitch)*cos(roll) + Zc*cx*cos(pitch)*cos(yaw) - Xc*cy*cos(pitch)*sin(roll) - Zc*u*cos(pitch)*cos(yaw) - Zc*cy*cos(roll)*sin(yaw) + Xc*v*cos(pitch)*sin(roll) - cx*z0*cos(pitch)*cos(yaw) - Zc*f*sin(roll)*sin(yaw) + Zc*v*cos(roll)*sin(yaw) + u*z0*cos(pitch)*cos(yaw) + cy*z0*cos(roll)*sin(yaw) - Zc*f*cos(roll)*cos(yaw)*sin(pitch) + Zc*cy*cos(yaw)*sin(pitch)*sin(roll) + f*z0*cos(roll)*cos(yaw)*sin(pitch) - Zc*v*cos(yaw)*sin(pitch)*sin(roll) - cy*z0*cos(yaw)*sin(pitch)*sin(roll) + v*z0*cos(yaw)*sin(pitch)*sin(roll))./(cx*sin(pitch) - u*sin(pitch) + f*cos(pitch)*cos(roll) - cy*cos(pitch)*sin(roll) + v*cos(pitch)*sin(roll));
Yw=(Yc*cx*sin(pitch) - Yc*u*sin(pitch) - f*z0*cos(yaw)*sin(roll) + v*z0*cos(roll)*cos(yaw) + u*z0*cos(pitch)*sin(yaw) + Yc*f*cos(pitch)*cos(roll) + Zc*cy*cos(roll)*cos(yaw) - Yc*cy*cos(pitch)*sin(roll) + Zc*cx*cos(pitch)*sin(yaw) + Zc*f*cos(yaw)*sin(roll) - Zc*v*cos(roll)*cos(yaw) + Yc*v*cos(pitch)*sin(roll) - Zc*u*cos(pitch)*sin(yaw) - cy*z0*cos(roll)*cos(yaw) - cx*z0*cos(pitch)*sin(yaw) + v*z0*sin(pitch)*sin(roll)*sin(yaw) - Zc*f*cos(roll)*sin(pitch)*sin(yaw) + Zc*cy*sin(pitch)*sin(roll)*sin(yaw) + f*z0*cos(roll)*sin(pitch)*sin(yaw) - Zc*v*sin(pitch)*sin(roll)*sin(yaw) - cy*z0*sin(pitch)*sin(roll)*sin(yaw))./(cx*sin(pitch) - u*sin(pitch) + f*cos(pitch)*cos(roll) - cy*cos(pitch)*sin(roll) + v*cos(pitch)*sin(roll));
Zw=z0*ones(size(u));
s=-(f*(Zc - z0))./(cx*sin(pitch) - u*sin(pitch) + f*cos(pitch)*cos(roll) - cy*cos(pitch)*sin(roll) + v*cos(pitch)*sin(roll));

y = [Xw Yw Zw s]';
y = y(:);

end
