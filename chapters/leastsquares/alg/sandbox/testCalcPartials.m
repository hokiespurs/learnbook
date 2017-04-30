h = eps^(1/3);
x = [1 2 3;7 8 9]';
b = [.1 0.2]';

modelfun = @(b,x)(b(1)*sin(x(:,1))+x(:,2)*b(2));
derivfunB = @(b,x)([b(1)*cos(x(:,1)) b(2)*ones(size(x(:,1)))]);
derivfunJ = @(b,x)([sin(x(:,1)) x(:,2)]);

yhat = modelfun(b,x);

Bfun = @(xn)(modelfun(b,xn));
Jfun = @(xn)(modelfun(xn,x));

fprintf('BFUN\n');
dfdx = calcPartials(Bfun,x,h)
derivfunB(b,x)-dfdx

fprintf('JFUN\n');
dfdb = calcPartials(Jfun,b',h)
derivfunJ(b,x)-dfdb

%% Linear
h = eps^(1/3);
x = [1 2 3;7 8 9]';
b = [.1 0.5]';

modelfun = @(b,x)(b(1)*(x(:,1))+x(:,2)*b(2));
derivfunB = @(b,x)([b(1)*(x(:,1)) b(2)*ones(size(x(:,1)))]);
derivfunJ = @(b,x)([(x(:,1)) x(:,2)]);

yhat = modelfun(b,x);

Bfun = @(xn)(modelfun(b,xn));
Jfun = @(xn)(modelfun(xn,x));

fprintf('BFUN\n');
dfdxn = calcPartials(Bfun,x,h)
derivfunB(b,x)

fprintf('JFUN\n');
dfdxn = calcPartials(Jfun,b',h)
derivfunJ(b,x)

%% Do linear different
