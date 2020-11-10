function [U,V,X,C,S] = newgsvd(A,B,arg3)
%   GSVD   Generalized Singular Value Decomposition.
%   [U,V,X,C,S] = NEWGSVD(A,B) returns unitary matrices U and V,
%   a (usually) square matrix X, and nonnegative diagonal matrices
%   C and S so that
%
%       A = U*C*X'
%       B = V*S*X'
%       S'*S = I 
%
%   A and B must have the same number of columns, but may have
%   different numbers of rows.  If A is m-by-p and B is n-by-p, then
%   U is m-by-m, V is n-by-n and X is p-by-q where q = min(m+n,p).
%
%   SIGMA = GSVD(A,B) returns the vector of generalized singular
%   values, sqrt(diag(C'*C)./diag(S'*S)).
%
%   The nonzero elements of S are always on its main diagonal.  If
%   m >= p the nonzero elements of C are also on its main diagonal.
%   But if m < p, the nonzero diagonal of C is diag(C,p-m).  This
%   allows the diagonal elements to be ordered so that the generalized
%   singular values are nondecreasing.
%   
%   GSVD(A,B,0), with three input arguments and either m or n >= p,
%   produces the "economy-sized" decomposition where the resulting
%   U and V have at most p columns, and C and S have at most p rows.
%   The generalized singular values are diag(C)./diag(S).
%   
%   When I = eye(size(A)), the generalized singular values, gsvd(A,I),
%   are equal to the ordinary singular values, svd(A), but they are
%   sorted in the opposite order.  Their reciprocals are gsvd(I,A).
%
%   In this formulation of the GSVD, no assumptions are made about the
%   individual ranks of A or B.  The matrix X has full rank if and only
%   if the matrix [A; B] has full rank.  In fact, svd(X) and cond(X) are
%   are equal to svd([A; B]) and cond([A; B]).  Other formulations, eg.
%   G. Golub and C. Van Loan, "Matrix Computations", require that null(A)
%   and null(B) do not overlap and replace X by inv(X) or inv(X').
%   Note, however, that when null(A) and null(B) do overlap, the nonzero
%   elements of C and S are not uniquely determined.
%
%   Class support for inputs A,B:
%      float: double, single
%
%   See also SVD.
