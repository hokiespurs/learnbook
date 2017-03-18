%% Example Nonlinear Least Squares Script
%% read data
data=csvread('ts.csv');
t = data(:,1);                  % raw time observations
y = data(:,2);                  % raw elevation observations
stdy = data(:,3);               % reported std for each observation
W = inv(diag(stdy.^2));         % Create a Weight Matrix based on the stdy

%% Solve Least Squares
X = [1.5; 1];                   % initial unknowns guess

So2 = inf; dSo2 = 1; iter = 0;  % initialize for while loop

m = NPTS;                       % number of observations
n = 2;                          % number of unknowns
dof = m-n;                      % degrees of freedom
while dSo2>=0 && iter<100 %loop until So2 increases or exceed 100 iteration
    J = [sin(2*pi/2.*t + X(2)) X(1)*cos(2*pi/2.*t + X(2))];
    K = y - X(1)*sin(2*pi/2.*t + X(2));
    dX = (J'*W*J)\J'*W*K;       % Loop Delta Estimate
    X = X + dX;                 % Loop Estimate
    V = J*dX-K;                 % Residuals
    dSo2 = So2 - V'*W*V/dof;    % Change in Reference Variance
    So2 = (V'*W*V)/dof;         % Reference Variance
    iter = iter + 1;
end
Q = inv(J'*W*J);                % cofactor
Sx = So2 * Q;                   % covariance of unknowns
Sl = J * Sx * J';               % covariance of observations
stdX = sqrt(diag(Sx));          % std of solved unknowns
Lhat = J * X;                   % predicted L values
RMSE = sqrt(V'*V/m);            % RMSE

% Perform chi^2 Test
alphaval = 0.05;
chi2 = dof * So2;
chi2low = chi2inv(alphaval/2,dof);
chi2high = chi2inv(1-alphaval/2,dof);
% Recalculate Values with So2 = 1 because CHI2 passed
So2_b = 1;                      % set reference variance equal to 1
Sx_b = So2_b * Q;               % covariance of unknowns     w/ So2 = 1
Sl_b = J * Sx_b * J';           % covariance of observations w/ So2 = 1
stdX_b = sqrt(diag(Sx_b));      % std of solved unknowns     w/ So2 = 1

%% Using Matlab Built in Function NLINFIT
beta0 = [1.5; 1];
modelfun = @(b,x)(b(1)*sin(2*pi/2*x+b(2)));
[mat_X,mat_V,mat_J,mat_Sx,mat_So2,ErrorModelInfo] = ...
    nlinfit(t,y,modelfun,beta0,'Weights',1./noisescale.^2);

%% Print chi2 tests
if chi2>chi2low || chi2<chi2high
    fprintf('Chi2 Test Passed, We cant say So2 doesnt equal 1\n');
else
    fprintf('Chi2 Test Failed, So2 not equal to 1\n');
end

texstr = printLatexGoodnessOfFit(So2,dof,alphaval);
fid = fopen('../tab/exampleNonlinearLSresultschi2.tex','w+t');
fprintf(fid,texstr);
fclose(fid);

%% Plot Raw Data
f = figure(1);clf;
f.Position = [415 346 1110 384];
plot(t,y,'b.','markersize',10)
axis equal
xlim([0 MAXT]);
grid on
xlabel('time(s)','fontsize',20,'interpreter','latex')
ylabel('elevation(m)','fontsize',20,'interpreter','latex')
title('Raw Timeseries Observations','fontsize',24,'interpreter','latex')
saveas(f,'../fig/nonlineardata.png')

%% Plot Results
hold on
A = X(1);
PHI = X(2);
xi = linspace(0,MAXT,10000);
yi = A * sin(2*pi/2*xi + PHI);
plot(xi,yi,'g-')
legend({'Raw Observations','Nonlinear Least Squares'},...
    'fontsize',16,'interpreter','latex');
title('Nonlinear Least Squares','fontsize',24,'interpreter','latex')
saveas(f,'../fig/nonlinearresult.png')

%% Make Latex Results Table
fid = fopen('../tab/exampleNonlinearLSresultsTable.tex','w+t');
T = cell(3,3);

T{1,1} = sprintf('$n = %.0f$',n);
T{1,2} = sprintf('$m = %.0f$',m);
T{1,3} = sprintf('$dof = %.0f$',dof);

T{2,1} = ['$\\hat{X} = $' printLatexMatrix(X,'%.2f')];

T{2,2} = sprintf('$S_0^2 = %.2f$ ',So2);

T{2,3} = ['$Q_{xx} = $ ' printLatexMatrix(Q,'%0.2f')];

T{3,1} = ['$\\Sigma_{xx} = $ ' printLatexMatrix(Sx,'%0.4f')];

T{3,2} = ['$\\sigma_{\\hat{X}} = $ ' printLatexMatrix(stdX,'%0.4f')];

T{3,3} = sprintf('$RMSE = %.2f$',RMSE);

fprintf(fid,printLatexTable(T));
fclose(fid);

%% Make Latex Results Table after CHI2 test
fid = fopen('../tab/exampleNonlinearLSresultsTable2.tex','w+t');

T2 = cell(1,3);

T2{1,1} = sprintf('$S_0^2 = %.2f$ ',So2_b);
T2{1,2} = ['$\\Sigma_{xx} = $ ' printLatexMatrix(Sx_b,'%0.4f')];
T2{1,3} = ['$\\sigma_{\\hat{X}} = $ ' printLatexMatrix(stdX_b,'%0.4f')];

fprintf(fid,printLatexTable(T2));
fclose(fid);