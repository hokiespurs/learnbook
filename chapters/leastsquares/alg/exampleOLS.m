%% Example Ordinary Least Squares
%% observations
x = [0 1 2 3 4];
y = [0 1 7 13 24];
%% Define A and L based on observation equation
A = [x(:).^2 x(:) ones(size(x(:)))];
L = y(:);
%% Calculate Ordinary Least Squares
% X = inv(A'*A)*A'*L; %Replaced 'inv(A)*b' with 'A\b' per Matlab docs
X = (A'*A)\A'*L;                % unknowns
V = A * X - L;                  % residuals
dof = numel(L) - size(A,2);     % degrees of freedom
So2 = V'*V/dof;                 % Reference Variance
Q = inv(A'*A);                  % cofactor
Sx = So2.^2 * Q;                % covariance
Lhat = A * X;                   % predicted L values
r2 = var(Lhat)/var(L);          % R^2 Skill

%% Plot Results
f = figure(1);clf
xx = -1:0.01:5;
yy = X(1)*xx.^2+X(2)*xx+X(3);
plot(xx,yy,'k','linewidth',3);
hold on
plot(x,y,'r.','markersize',30);
grid on
xlabel('X','fontsize',20,'interpreter','latex')
ylabel('Y','fontsize',20,'interpreter','latex')
title('Ordinary Least Squares','fontsize',24,'interpreter','latex')
h = text(-0.5,25,'y = ax$^2$ + bx + c','interpreter','latex');
h.HorizontalAlignment = 'left';
h.FontSize = 20;
h.VerticalAlignment = 'bottom';
h.BackgroundColor = 'w';
astr = sprintf('a = %.2f $\\pm$ %.2f',X(1),sqrt(Sx(1,1)));
bstr = sprintf('b = %.2f $\\pm$ %.2f',X(2),sqrt(Sx(2,2)));
cstr = sprintf('c = %.2f $\\pm$ %.2f',X(3),sqrt(Sx(3,3)));
h = text(0,25,{astr,bstr,cstr},'interpreter','latex');
h.FontSize = 16;
h.BackgroundColor = 'w';
h.VerticalAlignment = 'top';
hl=legend({'Ordinary Least Squares Fit','Observation Data'},'fontsize',16);
hl.Interpreter = 'latex';
hl.Location = 'northwest';
saveas(f,'OLSexample.png');