%% CONSTANTS for plane fitting Nx(x-x0)+Ny(y-y0)+Nz(z-z0) = 0
rng(10);
NOBS = 50;
XRANGE = [0 100];
YRANGE = [0 100];
TRUEB = [2 1 5]';

NOISEMAGNITUDEX = 10;
NOISEMAGNITUDEY = 10;
NOISEMAGNITUDEZ = 20;

XNOISEREPORTEDSTDNOISE = 0.1;
YNOISEREPORTEDSTDNOISE = 0.1;
ZNOISEREPORTEDSTDNOISE = 0.1;

MINCOMPUTEDSTD = 0.01;
GUESSB = [1 1 7]';
GUESSBCOV = diag([0.1 1]);

PLOTX = [0 0 100 100];
PLOTY = [0 100 100 0];
%% Compute Plane Data
% beta = [a b c d]  % ax + by + c = z
planeval0 = @(b,x)(b(1)*x(:,1) + b(2)*x(:,2) + b(3) - x(:,3));
planezval = @(b,x)(b(1)*x(:,1) + b(2)*x(:,2) + b(3));

xtrue = rand(NOBS,1)*(XRANGE(2)-XRANGE(1)) + XRANGE(1);
ytrue = rand(NOBS,1)*(YRANGE(2)-YRANGE(1)) + YRANGE(1);
ztrue = planezval(TRUEB,[xtrue ytrue]);

xerrorest = abs(normrnd(0,NOISEMAGNITUDEX,NOBS,1));
yerrorest = abs(normrnd(0,NOISEMAGNITUDEY,NOBS,1));
zerrorest = abs(normrnd(0,NOISEMAGNITUDEZ,NOBS,1));

xerrorest(xerrorest<MINCOMPUTEDSTD) = MINCOMPUTEDSTD;
yerrorest(yerrorest<MINCOMPUTEDSTD) = MINCOMPUTEDSTD;
zerrorest(zerrorest<MINCOMPUTEDSTD) = MINCOMPUTEDSTD;

xnoise = normrnd(0,xerrorest,NOBS,1);
ynoise = normrnd(0,yerrorest,NOBS,1);
znoise = normrnd(0,zerrorest,NOBS,1);

x = xtrue + xnoise;
y = ytrue + ynoise;
z = ztrue + znoise;

%% Linear unweighted
X = [x y];
Y = z;
[betacoef,R,J,CovB,MSE,ErrorModelInfo] = lsr2(X,Y,planezval,GUESSB,'verbose',true);

%plot
f5=figure(5);clf;hold on
f5.Name = 'Nonlinear Unweighted';
set(f5,'WindowStyle','docked')
% for i=1:numel(x)
%     plot3([x(i) xtrue(i)],[y(i) ytrue(i)],[z(i) ztrue(i)],'b.-')
% end
plot3(x,y,z,'b.')

PLOTZ = planezval(TRUEB,[PLOTX' PLOTY']);
fill3(PLOTX,PLOTY,PLOTZ,'g')
alpha 0.2
CALCZ = planezval(betacoef,[PLOTX' PLOTY']);
fill3(PLOTX,PLOTY,CALCZ,'r')
alpha 0.2

%% Nonlinear unweighted
X = [x y z];
Y = zeros(size(x));
[betacoef,R,J,CovB,MSE,ErrorModelInfo] = lsr2(X,Y,planeval0,GUESSB,'type','nonlinear','verbose',true);

%plot
f5=figure(5);clf;hold on
f5.Name = 'Nonlinear Unweighted';
set(f5,'WindowStyle','docked')
% for i=1:numel(x)
%     plot3([x(i) xtrue(i)],[y(i) ytrue(i)],[z(i) ztrue(i)],'b.-')
% end
plot3(x,y,z,'b.')

PLOTZ = planezval(TRUEB,[PLOTX' PLOTY']);
fill3(PLOTX,PLOTY,PLOTZ,'g')
alpha 0.2
CALCZ = planezval(betacoef,[PLOTX' PLOTY']);
fill3(PLOTX,PLOTY,CALCZ,'r')
alpha 0.2

%% Total
[betacoef,R,J,CovB,MSE,ErrorModelInfo] = lsr2(X,Y,planeval0,GUESSB,...
    'type','total','verbose',true);
% test plot
CALCZ = planezval(betacoef,[PLOTX' PLOTY']);
fill3(PLOTX,PLOTY,CALCZ,'y')
alpha 0.2

