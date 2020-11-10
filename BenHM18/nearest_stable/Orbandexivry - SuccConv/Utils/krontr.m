function [M] = krontr(X,Y)
% KRONTR(X,Y) is the Kronecker tensor product of X and Y where the columns
% are permuted. The permutation is done so that
%        krontr(X,Y)*vec(A)= kron(X,Y)*vec(A.')
%
% see KRON for details on the usual Kronecker product
%
%Copyright 2009: François-Xavier Orban de Xivry for
% Universite Catholique de Louvain
% Avenue Georges Lemaitre,4-6
% 1348 Louvain-la-Neuve
% Belgium

M = kron(X,Y);
E = permut(size(M,1));

M = M*E;% E==E'et E*E' = I donc M*E == M*E' 


end

%%

function [E]= permut(n)
%PERMUT(n) produces a (unitary) permutation matrix E of size n by n.
%
%When premultiplied (resp. postmultiplied) the permutation is built so that
%E*vec(A.') = vec(A) (resp. vec(A.').'*E = vec(A).') for any matrix A of
%size sqrt(n) by sqrt(n), i.e it permutes the lines (resp.columns) of any
%matrix. Thus, M*vec(A.') = (M*E')*vec(A)
%
%Input : n : the dimension of the square matrix E
%output: E : The permutation matrix
%
%Copyright 2009: François-Xavier Orban de Xivry for
% Universite Catholique de Louvain
% Avenue Georges Lemaitre,4-6
% 1348 Louvain-la-Neuve
% Belgium

if n <= 0
    error('The dimension of the permutation matrix is negative')
end

%dimension of the matrix A in its square form
dim = sqrt(n);

%initializing the permutation matrix
E = speye(n);

for i = 1 : dim
    for j = 1 : i

        if i ~= j
            %index of the lines that needs to be permuted
            k = (i-1)*dim + j;
            l = (j-1)*dim + i;

            E(k,k)= 0;
            E(k,l) = 1;
            E(l,k) = 1;
            E(l,l) = 0;

        end
    end
end
end