function texstr = printLatexCovariance(Xstr)
%%

N = numel(Xstr);

X = cell(N,N);

for i=1:N
    for j=1:N
        if i==j
           X{i,j}=['\\sigma_{' Xstr{i} '}^2'];
        else
           X{i,j}=['\\sigma_{' Xstr{i} Xstr{j} '}'];
        end
    end
end

texstr = '\\[\n\\Sigma = \n\\begin{bmatrix}\n';
for i=1:N
    for j=1:N
            texstr = [texstr X{i,j}];
        if j~=N
           texstr = [texstr sprintf(' & ')];            
        end
    end
    texstr = [texstr ' \\\\\n'];
end
texstr = [texstr '\\end{bmatrix}\n'];
texstr = [texstr '\\]\n'];

fprintf(texstr)
