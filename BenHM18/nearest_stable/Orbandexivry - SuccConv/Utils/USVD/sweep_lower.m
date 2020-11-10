function [A,B,U,V,Q] = sweep_lower(A,B)
%performs a sweep of the usvd either 'lower' or 'upper' according to the
%fact that A' and B are both square upper or lower  triangular matrices.


U = eye(n);
V = eye(n);
Q = eye(n);

for i = 1:n-1
    for j = i+1:n
        %computing the svd of a small 2 by 2 matrix
        E = A(:,[i j])'*B(:,[i j]);
        [d,u,v] = diag2by2(E);
        
        %applying rotation to each term of the product
        A(:,[i j]) = A(:,[i j])*u;
        B(:,[i j]) = B(:,[i j])*v;
        
        %killing the upper off diagonal term that remains in B (and A)
        [q] = kill_offdiag(B([i j],[i j]),'lower');
        
        A([i j],:) = q*A([i j],:);
        B([i j],:) = q*B([i j],:);
        
        %computing overall rotation
        U(:,[i j]) = U(:,[i j ])*u;
        V(:,[i j]) = V(:,[i j])*v;
        Q([i j],:) = q*Q([i j],:);
    end
end