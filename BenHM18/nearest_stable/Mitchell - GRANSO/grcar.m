% grcar(n,k) returns an n-by-n Toeplitz matrix with -1s on the subdiagonal, 
% 1s on the diagonal, and k superdiagonals of 1s. The default is k = 3. 

function A = grcar(n,k) 

A = eye(n); 
A = A - diag(ones(n-1,1),-1); 

for i = 1 : k
    A = A + diag(ones(n-i,1),i); 
end
