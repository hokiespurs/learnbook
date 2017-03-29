%% Example Ordinary Least Squares Script
%% observations
x = [0 1 2 3 4];
y = [5 1 7 13 24];

%% Define A and L based on observation equation
A = [x(:).^2 x(:) ones(size(x(:)))];
L = y(:);

%% Calculate Ordinary Least Squares
m = numel(L);                   % number of observations
n = size(A,2);                  % number of unknowns
dof = m-n;                      % degrees of freedom
X = (A'*A)\A'*L;                % unknowns(inv(A)*b' with 'A\b' per Matlab)
V = A * X - L;                  % residuals
So2 = V'*V/dof;                 % Reference Variance
Q = inv(A'*A);                  % cofactor
Sx = So2 * Q;                   % covariance of unknowns
Sl = A * Sx * A';               % covariance of observations
stdX = sqrt(diag(Sx));          % std of solved unknowns
Lhat = A * X;                   % predicted L values
r2 = var(Lhat)/var(L);          % R^2 Skill
RMSE = sqrt(V'*V/m);            % RMSE

%% Using Matlab Built in Function LSCOV
[matlab_X, matlab_stdX, matlab_So2, matlab_Sx] = lscov(A,L); %same results

%% Example Using LSRLIN
[X,Sx,lsainfo] = lsrlin(A,L);
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
saveas(f,'../fig/OLSexample.png');

%% Make Latex Results Table
fid = fopen('../tab/exampleOLSresultsTable.tex','w+t');
T = cell(4,3);

T{1,1} = sprintf('$n = %.0f$',n);
T{1,2} = sprintf('$m = %.0f$',m);
T{1,3} = sprintf('$dof = %.0f$',dof);

T{2,1} = ['$\\hat{X} = $' printLatexMatrix(X,'%.2f')];

T{2,2} = ['$V = $ ' printLatexMatrix(V,'%0.2f')];

T{2,3} = sprintf('$S_0^2 = %.2f$ ',So2);

T{4,1} = ['$Q_{xx} = $ ' printLatexMatrix(Q,'%0.2f')];

T{3,1} = ['$\\Sigma_{xx} = $ ' printLatexMatrix(Sx,'%0.2f')];

T{3,3} = ['$\\Sigma_{\\hat{l}\\hat{l}} = $ ' printLatexMatrix(Sl,'%0.2f')];

T{3,2} = ['$\\sigma_{\\hat{X}} = $ ' printLatexMatrix(stdX,'%0.2f')];

T{4,2} = ['$\\hat{L} = $' printLatexMatrix(Lhat,'%0.2f')];

T{4,3} = sprintf('$R^2 = %.4f$ \\\\hspace{1cm} $RMSE = %.2f$',r2,RMSE);

fprintf(fid,printLatexTable(T));
fclose(fid);
      