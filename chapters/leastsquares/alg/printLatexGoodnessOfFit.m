function texstr = printLatexGoodnessOfFit(So2,dof,significance)

chi2 = dof * So2;
chi2low = chi2inv(significance/2,dof);
chi2high = chi2inv(1-significance/2,dof);
if chi2>chi2low
   passstr1 = '(PASS)'; 
else
    passstr1 = '(FAIL)';
end
if chi2<chi2high
   passstr2 = '(PASS)'; 
else
    passstr2 = '(FAIL)';
end

if chi2>chi2low && chi2<chi2high
    resultstr = sprintf('The results pass the $\\\\chi^2$ Goodness of Fit Test at the %.2f significance level, so redo the final calculations with $S_0^2$ equal to 1.',significance);
else
    resultstr = sprintf('The results fail the $\\\\chi^2$ Goodness of Fit Test at the %.2f significance level, so check for blunders or a poor stochastic model',significance);
end

texstr = 'A two tailed $\\chi^2$ Goodness of Fit Test is performed at the';
texstr = [texstr sprintf(' %.2f',significance) ' Significance level:\n'];
texstr = [texstr '\\begin{align*}\n'];
texstr = [texstr '\\text{Null Hypothesis } H_0 &: S_0^2 = 1 \\\\\n'];
texstr = [texstr '\\text{Alternative Hypothesis } H_a &: S_0^2 \\neq 1 \\\\\n'];
texstr = [texstr '\\text{Significance } &: \\alpha = ' sprintf('%.2f \n',significance)];
texstr = [texstr '\\end{align*}\n'];
texstr = [texstr 'Test Statistic:\n'];
texstr = [texstr '\\[\n'];
texstr = [texstr '\\chi^2 = \\dfrac{vS_0^2}{\\sigma^2} = \\dfrac{dof\\times S_0^2}{1}'];
texstr = [texstr ' = dof\\times S_0^2 = ' sprintf('%.2f',chi2) '\n'];
texstr = [texstr '\\]\n'];
texstr = [texstr 'Rejection Region:\n'];
texstr = [texstr '\\begin{align*}\n'];
texstr = [texstr sprintf('%.2f',chi2) ' = \\chi^2 &> \\chi_{(\\alpha/2,dof)}^2 = ' sprintf('%.2f',chi2low) ' \\hspace{1cm} &' sprintf('%s\\\\\\\\\n',passstr1)];
texstr = [texstr sprintf('%.2f',chi2) ' = \\chi^2 &< \\chi_{(1-\\alpha/2,dof)}^2 = ' sprintf('%.2f',chi2high) ' \\hspace{1cm} &' sprintf('%s\n',passstr2)];

texstr = [texstr '\\end{align*}\n'];
texstr = [texstr 'Results:\n\n'];
texstr = [texstr resultstr '\n'];

end