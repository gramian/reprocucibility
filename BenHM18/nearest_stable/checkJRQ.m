% This code allwos to solve 
% min_{ P PD, J=-J', R PSD } ||A P - (R-J)||
% using CVX available from http://cvxr.com/cvx/ 
% 
% 
% ****** Input ******
% A       : a square matrix  
% ****** Output ******
% (J,R,P) : P PD, J=-J', R PSD such that || A P - (J-R) || is minimized 
%  Q      : Q = inv(P) 
%  eq     : || A - (J-R) Q || 
%  ep     : || A P - (J-R) || 
% 
% ep = 0   IFF   the matrix A is stable. 

function [J,R,Q,P,eq,ep] = checkJRQ(A)

n = size(A,1); 

cvx_begin quiet sdp 
    variable P(n,n);
    variable J(n,n); 
    variable R(n,n); 
    
    minimize( norm( A*P - (J-R) , 'fro')  );
    subject to
        P - 1e-3*eye(n) == hermitian_semidefinite( n ); 
        R ==  hermitian_semidefinite( n ); 
        J == -J'; 
cvx_end
Q = inv(P); 
eq = norm( A - (J-R)*Q , 'fro'); 
ep = norm( A*P - (J-R) , 'fro'); 