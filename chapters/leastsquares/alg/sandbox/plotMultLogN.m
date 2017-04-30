function h = plotMultLogN(t,amp,mu,sigma,toffset,varargin)
%%
% t  : time vector
% amp: amplitudes
% mu : means
% sigma: sigmas
% 
% h : plot handles
%% Example
clf
t = 0:0.001:2.25;
amp = [3 1 10];
mu = [0.25 0.3 0];
sigma = [.06 0.05 .4];
toffset = [-.685 0.5 0];
% cmapBlue = [231 242 85;
%             40 227 53;
%             40 214 226]/255;
cmapBlue = [0 0.25 0.25;
            0 0.5 0.5;
            0 0.75 0.75];
cmapBlack = zeros(3);

FalloffT = 1.85;
SigmaT = .15;
legendNames = {'Surface Return','Bottom Return','Water Volume Backscatter'};
%% Optional parameters
documulative = true;
dofill = true;
cmapFace = parula(numel(amp));
cmapLine = parula(numel(amp));
cmapFace = cmapBlue;
cmapLine = cmapBlack;
plotHorizontally = true;
%% Constants
nWaveforms = numel(amp);
nTimesteps = numel(t);
%% Generate Individual Waveforms
w = nan(nTimesteps,nWaveforms);
for iWaveformNum = 1:nWaveforms
    iAmp = amp(iWaveformNum);
    iMu = mu(iWaveformNum);
    iSigma = sigma(iWaveformNum);
    iToffset = toffset(iWaveformNum);
    w(:,iWaveformNum) = iAmp * lognpdf(t-iToffset,iMu,iSigma);
end
%% HARD CODE FALLOFF LAST ONE
iW = w(:,end);
indT = t>=FalloffT;
multfact = ones(size(t));
multfact(indT)=normpdf(t(indT),FalloffT,SigmaT)./(max(normpdf(t(indT),FalloffT,SigmaT)));
iW = iW .* multfact';
w(:,end)=iW;
%% Cumulative Waveforms
if documulative
   w = fliplr(cumsum(fliplr(w),2));
end
%% Pad with 0s for plotting
wplot = [zeros(1,nWaveforms); w; zeros(1,nWaveforms)];
tplot = [t(1) t t(end)];

if plotHorizontally
%% Plot Waveforms Horizontally
h = nan(1,nWaveforms);
for iWaveformNum = 1:nWaveforms
    if dofill %make the plot area filled by a color
        hold on
        h(iWaveformNum) = area(tplot,wplot(:,iWaveformNum));
        set(h(iWaveformNum),'faceColor',cmapFace(iWaveformNum,:));
        set(h(iWaveformNum),'edgeColor',cmapLine(iWaveformNum,:));
    else
        
    end
end
% axes
xlabel('Time','interpreter','latex','fontsize',24)
ylabel('Amplitude','interpreter','latex','fontsize',24)
% title('Bathymetric Water Surface and Bottom Returns','interpreter','latex','fontsize',30)
else
%% Plot Waveforms Vertically
h = nan(1,nWaveforms);
for iWaveformNum = 1:nWaveforms
    if dofill %make the plot area filled by a color
        hold on
        h(iWaveformNum) = fill(wplot(:,iWaveformNum),tplot,cmapFace(iWaveformNum,:));
        set(h(iWaveformNum),'faceColor',cmapFace(iWaveformNum,:));
        set(h(iWaveformNum),'edgeColor',cmapLine(iWaveformNum,:));
    else
        
    end
end
% axes
set(gca,'ydir','reverse')
xlabel('Amplitude','interpreter','latex','fontsize',24)
ylabel('Time','interpreter','latex','fontsize',24)
% title('Bathymetric Water Surface and Bottom Returns','interpreter','latex','fontsize',30)
end
legend(legendNames,'fontsize',24,'interpreter','latex')
xticklabels('')
yticklabels('')
grid on
end