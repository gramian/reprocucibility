function [stock,Bmat,currentcol] = stockItLimited(pos,stock,Bmat,mat)

%computing the eigenvalues of mat
vp = eig(mat);

%we stock the value in stock
stock(pos,:) =  vp.';
Bmat(:,pos) = vec(mat);

if nargout == 3
currentcol = pos+1;
end


