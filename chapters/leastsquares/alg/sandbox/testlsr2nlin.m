function testlsr2nlin
%% Test Nonlinear Least Squares CONSTANTS
RANDSEED = 1;
TRUEBETA = [.2; .4];

NSAMPLES = 20;

XRANGE = [0 10];

NOISEMAGNITUDEX  = [0.05 0.5];
NOISEMAGNITUDEY  = [0.05 1];
NOISEMAGNITUDEXY = 10;

YEXPEQN = @(b,x) b(1)*exp(b(2)*x);
YEXPEQNTLS = @(b,x) b(1)*exp(b(2)*x(:,1)) - x(:,2);

PLOTELLIPSES = true;
PLOTERRORVECTORS = true;
OUTPUTVERBOSE = true;
EQUIVTHRESH = 1e-6;

%% Generate NonLinear Data
% y = ae^(bx)
% ae^(bx)-y = 0
rng(RANDSEED);
for i=1:NSAMPLES
    Sx(i)=rand(1)*diff(NOISEMAGNITUDEX) + NOISEMAGNITUDEX(1);
    Sy(i)=rand(1)*diff(NOISEMAGNITUDEY) + NOISEMAGNITUDEY(1);
    Sxy(i)=rand(1)*NOISEMAGNITUDEXY;
    S{i} = [Sx(i) Sxy(i);Sxy(i) Sy(i)];
    while ~isCovariance(S{i},2)
        Sxy(i)=randn(1)*NOISEMAGNITUDEXY;
        S{i} = [Sx(i) Sxy(i);Sxy(i) Sy(i)];
    end
    xyError = mvnrnd([0 0],S{i},1);
    xError(i) = xyError(1);
    yError(i) = xyError(2);
end
%    trueX = rand(NSAMPLES,1)*diff(XRANGE) + XRANGE(1);
trueX = linspace(XRANGE(1),XRANGE(2),NSAMPLES)';
trueY = YEXPEQN(TRUEBETA,trueX);
valX = trueX + xError';
valY = trueY + yError';

% Plot predicted Covariance Error Ellipse and Data Points
f = figure(2);clf;hold on
f.Name = 'Nonlinear';
set(f,'WindowStyle','docked')

ix = linspace(min(valX), max(valX),100);
plot(ix,YEXPEQN(TRUEBETA,ix),'k-');
if PLOTELLIPSES
    for i=1:NSAMPLES
        p = plotCovarianceFill(valX(i), valY(i),S{i});
        p.FaceColor = 'k';
        p.FaceAlpha = 0.2;
    end
end
if PLOTERRORVECTORS
    for i=1:NSAMPLES
        plot([trueX(i) valX(i)],[trueY(i) valY(i)],'r.-');
    end
end
plot(trueX,trueY,'b.','markersize',10)
plot(valX,valY,'rx','markersize',10)
axis equal
%% Test Nonlinear Least Squares
% Nonlinear Unweighted
[betacoef1,R,J,CovB,MSE,ErrorModelInfo] = lsr2(valX,valY,YEXPEQN,TRUEBETA,'verbose',OUTPUTVERBOSE);

%matlab check
[matlab_beta, R2, J2, CovB2, MSE2] = nlinfit(valX,valY,YEXPEQN,TRUEBETA);
% comparing values actually threw errors, lsr2 produces better RMSE
rmse_matlabway = sqrt(sum((YEXPEQN(matlab_beta,valX)-valY).^2))
rmse_lsrway = sqrt(sum((YEXPEQN(betacoef1,valX)-valY).^2))
rmse_matlabway - rmse_lsrway

% Nonlinear Weighted
w = 1./sqrt(Sx.^2+Sy.^2)';
[betacoef2,R,J,CovB,MSE,ErrorModelInfo] = lsr2(valX,valY,YEXPEQN,TRUEBETA,'verbose',OUTPUTVERBOSE,'type','nonlinear','weights',w);
%matlab check
[matlab_beta, R2, J2, CovB2, MSE2] = nlinfit(valX,valY,YEXPEQN,TRUEBETA,'Weights',w);
% matlab gives back weighted residuals, lsr2 doesnt, again lsr2 has better residuals

% Nonlinear Covariance
[betacoef3,R,J,CovB,MSE,ErrorModelInfo] = lsr2(valX,valY,YEXPEQN,TRUEBETA,'verbose',OUTPUTVERBOSE,'type','nonlinear','weights',diag(sqrt(Sx.^2+Sy.^2)));
%matlab check (MATLAB CANT INPUT COVARIANCE)
[matlab_beta, R2, J2, CovB2, MSE2] = nlinfit(valX,valY,YEXPEQN,TRUEBETA,'Weights',w);

% Total Least Squares
[betacoef4,R,J,CovB,MSE,ErrorModelInfo] = lsr2([valX valY],zeros(size(valY)),YEXPEQNTLS,TRUEBETA,'verbose',OUTPUTVERBOSE,'type','total','weights',blkdiag(S{:}));

% plots
cmap = jet(6);
p1 = plot(ix,YEXPEQN(betacoef1,ix),'-','color',cmap(1,:),'linewidth',3);
p2 = plot(ix,YEXPEQN(betacoef2,ix),'-','color',cmap(4,:),'linewidth',3);
p4 = plot(ix,YEXPEQN(betacoef4,ix),'-','color',cmap(6,:),'linewidth',3);

legend([p1 p2 p4],{'Unweighted','Weighted','TLS'},'fontsize',15,'location','best')

%% Test Explicit Jacobians YB and YX (Yields the same answer as TLS above)
Jyb = @(b,x) [exp(b(2)*x(:,1)) b(1).*x(:,1).*exp(b(2).*x(:,1))];
Jyx = @(b,x) bumphdiag([b(1)*b(2)*exp(b(2)*x(:,1)) -ones(size(x(:,1)))],1);
[betacoef4,R,J,CovB,MSE,ErrorModelInfo] = lsr2([valX valY],zeros(size(valY)),...
    YEXPEQNTLS,TRUEBETA,'verbose',OUTPUTVERBOSE,'type','total',...
    'weights',blkdiag(S{:}),'JacobianYB',Jyb,'JacobianYX',Jyx);

%% Test Max Iteration
[~,~,~,~,~,~] = lsr2(valX,valY,YEXPEQN,TRUEBETA,'verbose',OUTPUTVERBOSE,'maxiter',100);
[~,~,~,~,~,~] = lsr2(valX,valY,YEXPEQN,TRUEBETA,'verbose',OUTPUTVERBOSE,'maxiter',10);

%% Check DerivStep
[~,~,~,~,~,~] = lsr2(valX,valY,YEXPEQN,TRUEBETA,'verbose',OUTPUTVERBOSE,'derivstep',.001);

end

function compareVals(x,y,thresh)
if any(abs(x(:)-y(:))>thresh)
   error('values not equal'); 
end

end