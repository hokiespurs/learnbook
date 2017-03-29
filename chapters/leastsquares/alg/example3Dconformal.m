%example3DConformal
x = [1094.883 503.891 2349.343 1395.320];
y = [820.085 1598.698 207.658 1348.853];
z = [109.821 117.685 151.387 215.261];
X = [10037.810 10956.680 8780.080 10185.800];
Y = [5262.090 5128.170 4840.290 4700.210];
Z = [772.040 783.000 782.620 851.320];

sx = [7 11 6 5]*0.001;
sy = [8 8 5 8]*0.001;
sz = [5 9 7 9]*0.001;
sX = [5 4 2 3]*0.01;
sY = [6 6 4 5]*0.01;
sZ = [5 9 2 3]*0.01;

w = zeros(numel(x)*6,1);
w(1:6:end)=sx;
w(2:6:end)=sy;
w(3:6:end)=sz;
w(4:6:end)=sX;
w(5:6:end)=sY;
w(6:6:end)=sZ;

s = diag(w.^2); %square because std -> var for covariance matrix

Xo = calcXoEst(x,y,z,X,Y,Z); %calculate initial guess

Jfun = @(Xi)conformal3DJfun(Xi(1),Xi(2),Xi(3),Xi(4),Xi(5),Xi(6),Xi(7),...
                            x,y,z,X,Y,Z,logical([1 1 1 1 1 1 1]));
                        
Bfun = @(Xi)conformal3DBfun(Xi(1),Xi(2),Xi(3),Xi(4),Xi(5),Xi(6),Xi(7),...
                            x,y,z,X,Y,Z);
                        
Kfun = @(Xi)conformal3DKfun(Xi(1),Xi(2),Xi(3),Xi(4),Xi(5),Xi(6),Xi(7),...
                            x,y,z,X,Y,Z);

[X1,Sx1,lsainfo1] = lsrtls(Jfun,Kfun,Bfun,Xo,s);
%% example with scale locked to 0.95
myS = 0.9499549354;
Xo(1)=[];
Jfun = @(Xi)conformal3DJfun(myS,Xi(1),Xi(2),Xi(3),Xi(4),Xi(5),Xi(6),...
                            x,y,z,X,Y,Z,logical([0 1 1 1 1 1 1]));
                        
Bfun = @(Xi)conformal3DBfun(myS,Xi(1),Xi(2),Xi(3),Xi(4),Xi(5),Xi(6),...
                            x,y,z,X,Y,Z);
                        
Kfun = @(Xi)conformal3DKfun(myS,Xi(1),Xi(2),Xi(3),Xi(4),Xi(5),Xi(6),...
                            x,y,z,X,Y,Z);

[X2,Sx2,lsainfo2] = lsrtls(Jfun,Kfun,Bfun,Xo,s);

%%

function J = conformal3DJfun(S,rx,ry,rz,tx,ty,tz,x,y,z,X,Y,Z,flag)
m = numel(x);
J = nan(m*3,7);
J(1:3:end,1) = x.*cos(ry).*cos(rz) - z.*sin(ry) + y.*cos(ry).*sin(rz);
J(1:3:end,2) = 0;
J(1:3:end,3) = -S.*(z.*cos(ry) + x.*cos(rz).*sin(ry) + y.*sin(ry).*sin(rz));
J(1:3:end,4) = S.*cos(ry).*(y.*cos(rz) - x.*sin(rz));
J(1:3:end,5) = 1;
J(1:3:end,6) = 0;
J(1:3:end,7) = 0;
J(2:3:end,1) = y.*(cos(rx).*cos(rz) + sin(rx).*sin(ry).*sin(rz)) - x.*(cos(rx).*sin(rz) - cos(rz).*sin(rx).*sin(ry)) + z.*cos(ry).*sin(rx);
J(2:3:end,2) = S.*x.*(sin(rx).*sin(rz) + cos(rx).*cos(rz).*sin(ry)) - S.*y.*(cos(rz).*sin(rx) - cos(rx).*sin(ry).*sin(rz)) + S.*z.*cos(rx).*cos(ry);
J(2:3:end,3) = S.*sin(rx).*(x.*cos(ry).*cos(rz) - z.*sin(ry) + y.*cos(ry).*sin(rz));
J(2:3:end,4) = - S.*x.*(cos(rx).*cos(rz) + sin(rx).*sin(ry).*sin(rz)) - S.*y.*(cos(rx).*sin(rz) - cos(rz).*sin(rx).*sin(ry));
J(2:3:end,5) = 0;
J(2:3:end,6) = 1;
J(2:3:end,7) = 0;
J(3:3:end,1) = x.*(sin(rx).*sin(rz) + cos(rx).*cos(rz).*sin(ry)) - y.*(cos(rz).*sin(rx) - cos(rx).*sin(ry).*sin(rz)) + z.*cos(rx).*cos(ry);
J(3:3:end,2) = S.*x.*(cos(rx).*sin(rz) - cos(rz).*sin(rx).*sin(ry)) - S.*y.*(cos(rx).*cos(rz) + sin(rx).*sin(ry).*sin(rz)) - S.*z.*cos(ry).*sin(rx);
J(3:3:end,3) = S.*cos(rx).*(x.*cos(ry).*cos(rz) - z.*sin(ry) + y.*cos(ry).*sin(rz));
J(3:3:end,4) = S.*x.*(cos(rz).*sin(rx) - cos(rx).*sin(ry).*sin(rz)) + S.*y.*(sin(rx).*sin(rz) + cos(rx).*cos(rz).*sin(ry));
J(3:3:end,5) = 0;
J(3:3:end,6) = 0;
J(3:3:end,7) = 1;
J(:,~flag)=[]; % handle optional unknowns

end

function B = conformal3DBfun(S,rx,ry,rz,tx,ty,tz,x,y,z,X,Y,Z)
m = numel(x);

b = nan(3,6);
b(1:3:end,1) = S.*cos(ry).*cos(rz);
b(1:3:end,2) = S.*cos(ry).*sin(rz);
b(1:3:end,3) = -S.*sin(ry);
b(1:3:end,4) = -1;
b(1:3:end,5) = 0;
b(1:3:end,6) = 0;
b(2:3:end,1) = S.*cos(rz).*sin(rx).*sin(ry) - S.*cos(rx).*sin(rz);
b(2:3:end,2) = S.*(cos(rx).*cos(rz) + sin(rx).*sin(ry).*sin(rz));
b(2:3:end,3) = S.*cos(ry).*sin(rx);
b(2:3:end,4) = 0;
b(2:3:end,5) = -1;
b(2:3:end,6) = 0;
b(3:3:end,1) = S.*(sin(rx).*sin(rz) + cos(rx).*cos(rz).*sin(ry));
b(3:3:end,2) = S.*cos(rx).*sin(ry).*sin(rz) - S.*cos(rz).*sin(rx);
b(3:3:end,3) = S.*cos(rx).*cos(ry);
b(3:3:end,4) = 0;
b(3:3:end,5) = 0;
b(3:3:end,6) = -1;

B = kron(eye(m),b);
end

function K = conformal3DKfun(S,rx,ry,rz,tx,ty,tz,x,y,z,X,Y,Z)
m = numel(x);
K = nan(m*3,1);
K(1:3:end,1) = X - tx + S.*z.*sin(ry) - S.*x.*cos(ry).*cos(rz) - S.*y.*cos(ry).*sin(rz);
K(2:3:end,1) = Y - ty + S.*x.*(cos(rx).*sin(rz) - cos(rz).*sin(rx).*sin(ry)) - S.*y.*(cos(rx).*cos(rz) + sin(rx).*sin(ry).*sin(rz)) - S.*z.*cos(ry).*sin(rx);
K(3:3:end,1) = Z - tz - S.*x.*(sin(rx).*sin(rz) + cos(rx).*cos(rz).*sin(ry)) + S.*y.*(cos(rz).*sin(rx) - cos(rx).*sin(ry).*sin(rz)) - S.*z.*cos(rx).*cos(ry);
end

function Xo = calcXoEst(x,y,z,X,Y,Z)
%https://www.researchgate.net/publication/228608056_3-D_Coordinate_Transformations
rXYZ = sqrt((X-mean(X)).^2 + (Y-mean(Y)).^2 + (Z-mean(Z)).^2);
rxyz = sqrt((x-mean(x)).^2 + (y-mean(y)).^2 + (z-mean(z)).^2);
ScaleEst = std(rXYZ)/std(rxyz);

[R,rx,ry,rz] = calc3dRotEst(x,y,z,X,Y,Z);

rotxyz = R*([x(:) y(:) z(:)]');
Tshift = [mean(X) mean(Y) mean(Z)]-mean(rotxyz,2)';

Xo = [ScaleEst rx ry rz Tshift]';

end

function [R,rx,ry,rz] = calc3dRotEst(x,y,z,X,Y,Z)
%% estimate the 3d rotation matrix to go from xyz to XYZ
%https://www.researchgate.net/publication/228608056_3-D_Coordinate_Transformations
%% calculate mean

mX = mean(X);
mY = mean(Y);
mZ = mean(Z);

%% calc rotation estimate for each
R1 = calc3drot(x,y,z);
R2 = calc3drot(X,Y,Z);
R = R2'*R1;
%% decompose into thetas
rz = atan2(-R(2,1),R(1,1));
ry = acos(R(1,1)/cos(rz));
rx = acos(R(3,3)*cos(rz)/R(1,1));

end

function R = calc3drot(x,y,z)
mx = mean(x);
my = mean(y);
mz = mean(z);

a = [x(1) y(1) z(1)] - [mx my mz];
ah = a./norm(a);

b = [x(2) y(2) z(2)] - [mx my mz];
bh = b./norm(b);

theta = acosd(dot(ah,bh));

ph = cross(bh,ah)/sind(theta);

qh = cross(ah,ph);

R = [qh;ah;ph];

end