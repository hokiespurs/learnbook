function hx = bumphdiag(x,n)
hx = [];
for i=n:n:size(x,1)
    hx = blkdiag(hx,x(i-n+1:i,:));
end

end