%%
d = 100;
az = 45*pi/180;
sd = 10;
saz = 1*pi/180;

A = [cos(az) -d*sin(az); sin(az) d*cos(az)];
S = [10 0;0 1/180*pi].^2;
S2 = A*S*A';

x = d*cos(az);
y = d*sin(az);

f = figure(1);clf
hold on
plot([-100 100],[0 0],'-','color',[0.6 0.6 0.6])
plot([0 0],[-100 100],'-','color',[0.6 0.6 0.6])
plot([0 x],[0 y],'k-','linewidth',3);

plot(0,0,'b.','markersize',50)
plot(x,y,'g.','markersize',50);

text(10,5,'$\theta$','fontsize',40,'interpreter','latex')
text(40,50,'$d$','fontsize',40,'interpreter','latex')
text(75,70,'tree','fontsize',20)
text(1,-5,'observer','fontsize',20)

xlim([-10 100]);
ylim([-10 100]);
grid('minor')

xlabel('X (m)','fontsize',30,'interpreter','latex')
ylabel('Y (m)','fontsize',30,'interpreter','latex')

saveas(f,'../fig/treefig.png');

f = figure(1);clf
hold on
id = [90 110 110 90 90];
iaz = [46 46 44 44 46]*pi/180;
plot(id.*cos(iaz),id.*sin(iaz),'r-','linewidth',2)
p = plotCovarianceFill(x,y,enforceSymmetric(S2,1e-10));
p.FaceColor = 'w';
p.LineWidth = 2;
plot(x,y,'g.','markersize',50);
grid('minor')
xlabel('X (m)','fontsize',30,'interpreter','latex')
ylabel('Y (m)','fontsize',30,'interpreter','latex')
title('Tree Location Uncertainty','interpreter','latex','fontsize',36)
saveas(f,'../fig/treefigconfidence.png');

