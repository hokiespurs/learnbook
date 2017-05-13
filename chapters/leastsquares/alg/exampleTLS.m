%% Example TLS Script
%% Input data
x = [10 20 60 40 85];
y = [0 15 23 25 40];
S = blkdiag([45 -30;-30 30],[20 -10;-10 70],[80 4;4 4],[40 -13;-13 60],[30 -25;-25 30]);

%% Solve TLS
Xo =[0.5; 0]; 
X = Xo;                      % initial unknowns guess

So2 = inf; dSo2 = 1; iter = 0;     % initialize for while loop

m = numel(x);                      % number of observations
n = numel(X);                      % number of unknowns
dof = m-n;                         % degrees of freedom
while dSo2>0 && iter<100 %loop until So2 increases or exceed 100 iteration
    B = kron(eye(numel(y)),[X(1) -1]); %B
    J = [x(:) ones(size(x(:)))];   % Jacobian
    K = -(X(1)*x(:) + X(2) - y(:));% K Matrix
    Weq = inv(B*S*B');
    dX = (J'*Weq*J)\J'*Weq*K;      % Loop Delta Estimate
    X = X + dX;                    % Loop Estimate
    Veq = K;                       % Residuals
    dSo2 = So2 - Veq'*Weq*Veq/dof; % Change in Reference Variance
    So2 = (Veq'*Weq*Veq)/dof;      % Reference Variance
    iter = iter + 1;
end

V = S * B' * Weq * Veq;           % Observation Residuals
Q = inv(J'*Weq*J);                % cofactor
Sx = So2 * Q;                     % covariance of unknowns
Sl = J * Sx * J';                 % covariance of observations
stdX = sqrt(diag(Sx));            % std of solved unknowns
Lhat = J * X;                     % predicted L values
RMSE = sqrt(Veq'*Veq/m);          % RMSE

%% Example Using LSRTLS
Jfun = @(X)([x(:) ones(size(x(:)))]);
Kfun = @(X)(-(X(1)*x(:) + X(2) - y(:)));
Bfun = @(X)(kron(eye(numel(y)),[X(1) -1]));
[X2,Sx2,lsainfo] = lsrtls(Jfun,Kfun,Bfun,Xo,S);
%% Plot Results
xi = [0 100];
yiOLS = X(1)*xi + X(2);

f = figure(1);clf;
f.Position = [375         178        1000         667];
%plot Covariances
hold on

for i=1:numel(x)
    ind = (i-1)*2+1:(i-1)*2+2;
    plotCovarianceFill(x(i),y(i),S(ind,ind),0.5,inf,100,'faceColor',[0.2 0.2 0.2]);
    alpha 0.25
end
p3 = plot(-999,999,'o','MarkerFaceColor',[0.65 0.65 0.65],'MarkerEdgeColor','k','markersize',20);

p2 = plot(xi,yiOLS,'g-','linewidth',3);

%plot residuals
for i=1:numel(x)
   p4 = plot([x(i) x(i)],[y(i) y(i)-Veq(i)],'b-..','linewidth',3);
   p5 = plot([x(i) x(i)+V((i-1)*2+1)],[y(i) y(i)+V((i-1)*2+2)],'m-..','linewidth',3);
end
p1 = plot(x,y,'k.','markersize',30);


axis equal
ylim([-10 50]);
xlim([0 100])
xlabel('X','fontsize',20,'interpreter','latex')
ylabel('Y','fontsize',20,'interpreter','latex')
title('Total Least Squares','fontsize',24,'interpreter','latex')
hl=legend([p1 p2 p3 p4 p5],{'Raw Observations',...
    'Most Probable Line','Error Ellipse $\sigma_{50\%}$','Equivalent Residuals','Observation Residuals'},'fontsize',16);
hl.Interpreter = 'latex';
hl.Location = 'northwest';
grid on

% saveas(f,'../fig/TLSexampleA.png');
Xtls = X;
Vtls = V;
%% Make Latex Results Table
fid = fopen('../tab/exampleTLSresultsTable.tex','w+t');
T = cell(4,3);

T{1,1} = sprintf('$n = %.0f$',n);
T{1,2} = sprintf('$m = %.0f$',m);
T{1,3} = sprintf('$dof = %.0f$',dof);

T{2,1} = ['$\\hat{X} = $' printLatexMatrix(X,'%.2f')];

T{2,2} = ['$V_{eq} = $ ' printLatexMatrix(Veq,'%0.2f')];

T{2,3} = ['$V = $ ' printLatexMatrix(V,'%0.2f')];

T{4,1} = ['$Q_{xx} = $ ' printLatexMatrix(Q,'%0.2f')];

T{3,1} = ['$\\Sigma_{xx} = $ ' printLatexMatrix(Sx,'%0.2f')];

T{3,3} = ['$\\Sigma_{\\hat{l}\\hat{l}} = $ ' printLatexMatrix(Sl,'%0.2f')];

T{3,2} = ['$\\sigma_{\\hat{X}} = $ ' printLatexMatrix(stdX,'%0.2f')];

T{4,2} = ['$\\hat{L} = $' printLatexMatrix(Lhat,'%0.2f')];

T{4,3} = sprintf('$S_0^2 = %.2f \\\\hspace{1cm} RMSE = %.2f$',So2,RMSE);

fprintf(fid,printLatexTable(T));
fclose(fid);

%% Do General Least Squares to Compare

A = J;
L = y(:);
W = eye(5);

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

%% Plot Results
xi = [0 100];
yiOLS = X(1)*xi + X(2);
yiTLS = Xtls(1)*xi + Xtls(2);

f = figure(2);clf;
f.Position = [375         178        1000         667];
%plot Covariances
hold on
c1 = [228,26,28]/255;
c2 = [55,126,184]/255;
c3 = [228,26,28]/255;
c4 = [55,126,184]/255;

p2 = plot(xi,yiOLS,'-','color',c1,'linewidth',3);
p3 = plot(xi,yiTLS,'-','color',c2,'linewidth',3);

%plot residuals
for i=1:numel(x)
   p4 = plot([x(i) x(i)],[y(i) y(i)+V(i)],':','color',c3,'linewidth',2);
   p5 = plot([x(i) x(i)+Vtls((i-1)*2+1)],[y(i) y(i)+Vtls((i-1)*2+2)],':','color',c4,'linewidth',2);
end
p1 = plot(x,y,'k.','markersize',30);


axis equal
ylim([-10 50]);
xlim([0 100])
xlabel('X','fontsize',20,'interpreter','latex')
ylabel('Y','fontsize',20,'interpreter','latex')
title({'OLS vs TLS','$y = mx+b$'},'fontsize',24,'interpreter','latex')
hl=legend([p1 p2 p4 p3 p5],{'Raw Observations',...
    'OLS Fit','OLS Residual','TLS Fit','TLS Residual'},'fontsize',16);
hl.Interpreter = 'latex';
hl.Location = 'northwest';
grid on

saveas(f,'../fig/TLSexampleB.png');
