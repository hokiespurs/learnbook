%% Example General Least Squares Script
%% observations
xc = [1 2 3];
yc = [0 5 1];
x = [6 1 8];
y = [3 12 8];

Sc = [0.5 0.3 0 0 0 0;
      0.3 0.5 0 0 0 0;
      0 0 0.4 0.1 0 0;
      0 0 0.1 0.2 0 0;
      0 0 0 0 0.7 -0.4;
      0 0 0 0 -0.4 0.4]; %variance-covariance of control data
  
%% Define A, L, and W based on observation equation
A = nan(numel(x)*2,4);
A(1:2:end)=[x(:) -y(:) ones(size(x(:))) zeros(size(y(:)))];
A(2:2:end)=[y(:) x(:) zeros(size(y(:))) ones(size(x(:)))];

L = nan(numel(x)*2,1);
L(1:2:end) = xc;
L(2:2:end) = yc;

W = inv(Sc);

%% Calculate General Least Squares
m = numel(L);                   % number of observations
n = size(A,2);                  % number of unknowns
dof = m-n;                      % degrees of freedom
X = (A'*W*A)\A'*W*L;            % unknowns(inv(A)*b' with 'A\b' per Matlab)
V = A * X - L;                  % residuals
So2 = V'*W*V/dof;               % Reference Variance
Q = inv(A'*W*A);                % cofactor
Sx = So2 * Q;                   % covariance of unknowns
Sl = A * Sx * A';               % covariance of observations
stdX = sqrt(diag(Sx));          % std of solved unknowns
Lhat = A * X;                   % predicted L values
r2 = var(Lhat)/var(L);          % R^2 Skill
RMSE = sqrt(V'*V/m);            % RMSE

%% Using Matlab Built in Function LSCOV
[mat_X, mat_stdX, mat_So2, mat_Sx] = lscov(A,L,Sc); %same results

%% Plot Results
f = figure(1);clf
p1 = plot(xc,yc,'g.','markersize',20);
hold on
p2 = plot(x,y,'r.','markersize',20);
p3 = plot(Lhat(1:2:end),Lhat(2:2:end),'b.','markersize',20);
for i=1:3
    ind = (i-1)*2+1:(i-1)*2+2;
    plotCovarianceLine(xc(i),yc(i),Sc(ind,ind),0.5,inf,100,'g');
    p5 = plot([xc(i) x(i)],[yc(i) y(i)],'--','color',[0.6 0.6 0.6]);
end
p4 = plot(-999,999,'o','MarkerFaceColor','w','MarkerEdgeColor','g','markersize',20);
grid on
axis equal
xlim([-5 12]);
ylim([-1 13])
xlabel('X','fontsize',20,'interpreter','latex')
ylabel('Y','fontsize',20,'interpreter','latex')
title('GLS 2D Conformal Transformation','fontsize',24,'interpreter','latex')
hl=legend([p1 p2 p3 p4 p5],{'$(x_c,y_c)$',...
    '$(x,y)$','$(\hat{x_c},\hat{y_c})$',...
    'Control $\sigma_{50\%}$',...
    'correspondences'},'fontsize',16);
hl.Interpreter = 'latex';
hl.Location = 'northwest';
saveas(f,'../fig/GLSexample.png');

%% Make Latex Results Table
fid = fopen('../tab/exampleGLSresultsTable.tex','w+t');
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