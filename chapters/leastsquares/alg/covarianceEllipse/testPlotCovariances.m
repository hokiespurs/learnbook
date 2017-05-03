%% Test Covariance Plotting
%made up data
Sxx = 1;Syy = 2;Szz = 9;
Sxy = 1;Sxz = 1;Syz = 0.5;
C = [Sxx Sxy Sxz;Sxy Syy Syz;Sxz Syz Szz];

%MAG GPS
Cxx = 1.394500645924E-005;
Cyy = 2.371056244075E-005;
Czz = 9.695190117622E-005;
Cxy = 9.368919368385E-006;
Cxz = -2.005312258009E-005;
Cyz = -8.209277685364E-006;

C = [Cxx Cxy Cxz;
      Cxy Cyy Cyz;
      Cxz Cyz Czz];

CONFIDENCE = 0.8;DOF = 50;

f = figure(10);clf


% plotCovarianceLine + plotCovarianceRect2

subplot 221
plotCovarianceLine(0,0,C(1:2,1:2),CONFIDENCE,DOF,50,'linewidth',3);
axis equal
hold on
plotCovarianceRect2(0,0,C(1:2,1:2),CONFIDENCE,DOF,'k');
axis equal
xlabel('X','fontsize',20);
ylabel('Y','fontsize',20);

% plotCovarianceFill + plotCovarianceAx2
subplot 222
plotCovarianceFill(0,0,C(2:3,2:3),CONFIDENCE,DOF,50,'faceColor','c');
axis equal
hold on
plotCovarianceAxes2(0,0,C(2:3,2:3),CONFIDENCE,DOF,'linewidth',2);
axis equal
xlabel('Y','fontsize',20);
ylabel('Z','fontsize',20);

% plotCovarianceFill multiple levels
subplot 223
cmap = jet(5);
plotCovarianceFill(0,0,C([1 3],[1 3]),0.99,DOF,50,'faceColor',cmap(1,:));
hold on
plotCovarianceFill(0,0,C([1 3],[1 3]),0.95,DOF,50,'faceColor',cmap(2,:));
plotCovarianceFill(0,0,C([1 3],[1 3]),0.90,DOF,50,'faceColor',cmap(3,:));
plotCovarianceFill(0,0,C([1 3],[1 3]),0.75,DOF,50,'faceColor',cmap(4,:));
plotCovarianceFill(0,0,C([1 3],[1 3]),0.50,DOF,50,'faceColor',cmap(5,:));
plotCovarianceAxes2(0,0,C([1 3],[1 3]),0.99,DOF,'k','linewidth',2);
plotCovarianceRect2(0,0,C([1 3],[1 3]),0.99,DOF,'k');
axis equal
xlabel('X','fontsize',20);
ylabel('Z','fontsize',20);

% plotCovariance3surf + plotCovarianceRect3 + plotCovarianceAx3
subplot 224
plotCovarianceSurf(0,0,0,C,CONFIDENCE,DOF,[50 50],'faceColor','b');
alpha 0.1
hold on
plotCovarianceAxes3(0,0,0,C,CONFIDENCE,DOF,'linewidth',4);
plotCovarianceRect3(0,0,0,C,CONFIDENCE,DOF,'k');
axis equal
xlabel('X','fontsize',20);
ylabel('Y','fontsize',20);
zlabel('Z','fontsize',20);

%
figure(11);clf
problevels = [99.99 99.9 99 98 95 90 80 70 60 50 40 30 20 10 5 1];
cmap = parula(numel(problevels));
lstr = cell(1);
for i=1:numel(problevels)
    plotCovarianceFill(0,0,C([1 3],[1 3]),problevels(i)/100,DOF,50,...
        'faceColor',cmap(i,:));
    hold on
    lstr{i} = sprintf('%5.2f %%',problevels(i));
end
l = legend(lstr);
set(l,'fontsize',20)
plotRect([0 0],sqrt(diag(C([1 3],[1 3]))),'r','linewidth',3);
axis equal

xlabel('X','fontsize',20);
ylabel('Z','fontsize',20);