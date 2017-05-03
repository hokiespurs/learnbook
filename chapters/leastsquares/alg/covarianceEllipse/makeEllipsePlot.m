
C = [4 3; 3 16];
f = figure(1);clf;
f.Units = 'Normalize';
f.Position = [0.1 0.1 0.8 0.8];
problevels = [99.99 99.9 99 98 95 90 80 70 60 50 40 30 20 10 5 1];
cmap = parula(numel(problevels));
lstr = cell(1);
for i=1:numel(problevels)
    plotCovarianceFill(10,20,C,problevels(i)/100,inf,50,...
        'faceColor',cmap(i,:));
    hold on
    lstr{i} = sprintf('%5.2f %%',problevels(i));
end

% plotRect([0 0],sqrt(diag(C)),'r','linewidth',3);
axis equal
% lstr{i+1} = 'Standard Error Rectangle';
l = legend(lstr);
set(l,'fontsize',20)
l.Location = 'eastoutside';
xlim(10+[-10 10]);
ylim(20+[-20 20]);
xlabel('X','fontsize',50,'interpreter','latex');
ylabel('Y','fontsize',50,'interpreter','latex');
a
saveas(f,'../fig/ellipsefig.png');

f = figure(1);clf;
f.Units = 'Normalize';
f.Position = [0.1 0.1 0.8 0.8];


