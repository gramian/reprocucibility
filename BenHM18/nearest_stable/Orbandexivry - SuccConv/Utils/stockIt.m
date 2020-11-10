function [stock,Bmat] = stockIt(stock,Bmat,mat)

%dimensions
[nl,nc]= size(mat);

%if vector (polynomial) is in column transpose
if nl == 1
    mat = mat.';
    nl = nc;
    nc = 1;
end
%if polynomial creates a square matrix
if nc == 1
    Z = diag(ones(nl-1,1),-1);
    en = [zeros(nl-1,1);1];
    mat = Z - mat*en';
end

%computing the eigenvalues of mat
vp = eig(mat);

%we stock the value in stock
 stock = [stock; vp.'];
 Bmat = [Bmat vec(mat)];

