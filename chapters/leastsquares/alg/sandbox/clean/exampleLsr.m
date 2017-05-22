function exampleLsr



end


function [x,y,S,isoutlier]=generate2DSyntheticData(RANDSEED,TRUEBETA,...
    NSAMPLES,NOUTLIERS,OUTLIERMAG,OUTLIERSTD,XRANGE,NOISEMAGNITUDEX,...
    NOISEMAGNITUDEY,NOISEMAGNITUDEXY,MODELFUN)
%% Generate a random 2D dataset with covariances

% seed the random number generator
rng(RANDSEED);
%% Add Outliers
badX = rand(NOUTLIERS,1)*diff(XRANGE)+XRANGE(1);
badY = MODELFUN(TRUEBETA,badX) + (round(rand(NOUTLIERS,1))*2-1) * OUTLIERMAG + randn(NOUTLIERS,1)*OUTLIERSTD;

%% Generate NonLinear Data
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
trueY = MODELFUN(TRUEBETA,trueX);
valX = trueX + xError';
valY = trueY + yError';

% add outlines
x = [valX; badX];
y = [valY; badY];

isoutlier = true(size(x));
isoutlier(1:numel(valX))=false;
end