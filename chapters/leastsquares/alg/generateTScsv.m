%% Generate fake data for nonlinear least squares
TRUEAMP    = 1.425;
TRUEPERIOD = 2;
TRUEPHASE  = 0.95;
NPTS = 50;
MAXT = 10;
BASENOISE = 0.15;
SCALENOISE = 0.5;
RANDOMSEED = 1;

rng(RANDOMSEED);  %seed random number generator for consistent results
t = sort(rand(NPTS,1)*MAXT);
y = TRUEAMP * sin(2*pi/TRUEPERIOD .*t + TRUEPHASE);

%% add noise based on how far from 0 the data is
noisescale = abs(y) * SCALENOISE + BASENOISE;
% noisescale = ones(NPTS,1);
ynoise = noisescale.*randn(NPTS,1);
y = y + ynoise;
stdY = noisescale;

%% Save Data
csvwrite('ts.csv',[t y stdY]);