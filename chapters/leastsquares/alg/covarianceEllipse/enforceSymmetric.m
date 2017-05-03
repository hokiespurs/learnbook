function Y = enforceSymmetric(X,thresh)
Y=X;
for i=1:size(X,1)
   for j = 1:size(X,1)
      d = X(i,j)-X(j,i);
      if d>thresh
         error('Too Big of a difference to enforce symmetric'); 
      end
%       vals = [X(i,j) X(j,i)];
%       [~,lowind] = min(abs(vals));
%       Y(i,j)= vals(lowind);
      Y(j,i)= (X(i,j)+X(j,i))/2;
      Y(i,j)= Y(j,i);

   end
end

end