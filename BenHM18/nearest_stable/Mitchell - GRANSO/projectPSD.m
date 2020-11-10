% Project the matrix Q onto the PSD cone 
% 
% This requires an eigendecomposition and then setting the negative
% eigenvalues to zero. 

function Qp = projectPSD(Q) 

Q = (Q+Q')/2; 
if max(max(isnan(Q))) == 1 || max(max(isinf(Q))) == 1
    error('Input matrix has infinite of NaN entries');
end
[V,e] = eig(Q); 
Qp = V*diag(max(diag(e),0))*V'; 