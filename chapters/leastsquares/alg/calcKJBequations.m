function [J,B,K]=calcKJBequations(Eqs,invars,outvars)

clc;

K = Eqs;
% print anonymous function examples
fprintf('%%%% Anonymous Functions \n');
astr = ['(' sprintf('Xi(%.0f),',1:numel(outvars)) '...\n\t' ...
                 sprintf('%s,',invars)];
fprintf(['Jfun = @(Xi) Jfunction' astr '...\n\t' ...
    ' logical([' sprintf('%.0f ',ones(1,numel(outvars))) ...
    ']));\n']);

fprintf(['Bfun = @(Xi) Bfunction' astr(1:end-1) ');\n']);
fprintf(['Kfun = @(Xi) Kfunction' astr(1:end-1) ');\n']);

% print K, J, B matrices to screen
fprintf('K EQUATIONS \n');
fprintf('function K = Kfunction(');
tempstr = sprintf('%s,',[outvars invars]);
tempstr(end)=[];
fprintf([tempstr ')\n']);
fprintf('K = nan(numel(%s)*%.0f,1);\n',invars(1),numel(Eqs));
for iEq = 1:numel(Eqs)
   fprintf('K(%.0f:%.0f:end,1) = %s;\n',iEq,numel(Eqs),makeSymbolString(Eqs(iEq)));
end

fprintf('\nJ EQUATIONS \n');
fprintf('function J = Jfunction(');
fprintf('%s,',[outvars invars]);
fprintf('flag)\n');
fprintf('J = nan(numel(%s)*%.0f,%.0f);\n',invars(1),numel(Eqs),numel(outvars));
for iEq = 1:numel(Eqs)
    for jOut = 1:numel(outvars)
        J(iEq,jOut) = simplify(diff(-Eqs(iEq),outvars(jOut)));
        fprintf('J(%.0f:%.0f:end,%.0f) = %s;\n',iEq,numel(Eqs),jOut,makeSymbolString(J(iEq,jOut)));
    end
end
fprintf('J(:,~flag)=[]; %% handle optional unknowns\n');

fprintf('\nB EQUATIONS \n');
fprintf('function B = Bfunction(');
tempstr = sprintf('%s,',[outvars invars]);
tempstr(end)=[];
fprintf([tempstr ')\n']);
fprintf('B = nan(numel(%s)*%.0f,%.0f);\n',invars(1),numel(Eqs),numel(invars));
for iEq = 1:numel(Eqs)
    for jIn = 1:numel(invars)
        B(iEq,jIn) = simplify(diff(-Eqs(iEq),invars(jIn)));
        fprintf('B(%.0f:%.0f:end,%.0f) = %s;\n',iEq,numel(Eqs),jIn,makeSymbolString(B(iEq,jIn)));
    end
end
fprintf('B = bumphdiag(B,%.0f); %%this pads the matrix with 0s\n',numel(Eqs));

end

function strfix = makeSymbolString(S)
    strfix = sprintf('%s',S);
    strfix = strrep(strfix,'*','.*');
    strfix = strrep(strfix,'/','./');
    strfix = strrep(strfix,'^','.^');
end