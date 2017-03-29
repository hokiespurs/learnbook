x = -5:0.01:5;
n = [2 5 inf];
y = []; x95 = [];lstr=[];
for i=1:numel(n)
    y{i} = tpdf(x,n(i));
    x95{i} = tinv(0.05,n(i));
    lstr{i} = ['$t_{(\sigma=1,dof=' sprintf('%.0f',n(i)) ')}$' ];
end
% plot
f=figure(1);clf;hold on
set(gca,'fontsize',20)
cmap =  [        0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];

for i=1:numel(y)
    plot(x,y{i},'linewidth',3,'color',cmap(i,:))
end
for i=1:numel(y)
    plot([x95{i} x95{i} nan -x95{i} -x95{i}],[0 0.4 nan 0 0.4],'--','color',cmap(i,:),'linewidth',2)
    lstr{numel(y)+i} = sprintf('$95\\%% CI_{(dof=%.0f)} = %.2f\\times\\pm \\sigma$',n(i),-x95{i});
end

l = legend(lstr,'interpreter','latex','fontsize',24);
xlabel('X Value','fontsize',24,'interpreter','latex')
ylabel('Probability','fontsize',24,'interpreter','latex')
title('t-distribution with various degrees of freedom','fontsize',30,'interpreter','latex')

saveas(f,'../fig/tdist.png');