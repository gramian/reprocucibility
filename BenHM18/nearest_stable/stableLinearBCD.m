% Find the nearest stable matrix in the form (J-R)*Q to a given matrix A
% 
% ****** Description ******
% This code tries to solve the following optimization problem
% 
% min_{J,R,Q} || A - (J-R) * Q ||_F 
% 
% such that J = -J', R and Q PSD. 
%
% It uses a *block coordinate descent method* with block (J,R) and Q,
% otpimizing alternatively over (J,R) and Q. 
% 
% See the paper 
% On computing the distance to stability for matrices using linear
% dissipative Hamiltonian systems, Nicolas Gillis and Punit Sharma, 2016. 
% 
% ****** Input ******
% A        : an (unstable) matrix A
% maxiter  : the maximum number of iterations performed by the algorithm 
%            -default = 1e6. 
% timemax  : the maximum time alloted to the algorithm 
%            -default = 10. 
% (J,R,Q)  : initialization. 
%            -default: Q = I, J = (A-A')/2, R = projectPSD(-(A+A')/2). 
%
% ****** Output ******
% (J,R,Q)  : J=J', R and Q PSD, such that (J-R)*Q is close to A 
% e        : evolution of the approximation error (Frobenius norm)
% t        : cputime to compute e --plot(t,e) displays the error over time


function [J,R,Q,e,t] = stableLinearBCD(A,maxiter,timemax,J,R,Q); 

cput = cputime; 
n = size(A,1); 
if nargin <= 1
    maxiter = 1e6; 
end
if nargin <= 2
    timemax = 10; 
end
if nargin <= 3 || isempty(Q)
    Q = eye(n); 
    J = (A-A')/2; 
    R = projectPSD(J-A); 
    ejr =  norm( A - (J-R)*Q , 'fro'); 
else
   [J,R,ejr] = porthamJRfg(A,zeros(n),zeros(n),Q,10000); 
    ejr = ejr(end); 
end
% initial relative error
e(1) = ejr; 
t(1) = cputime - cput; 

% Scaling (J,R,Q) for better numerical stability 
[uq,Lq,vq] = svds(Q,1); 
Lq = Lq^2; 
% -- For Dq
JmR = J-R; 
JmRtJmR = JmR'*JmR; 
[ujmr,Ljmr,vjmr] = svds(JmRtJmR,1); 
% -- scaling such that Lq = Ljmr: a*Q -> a^2Lq et b b*JmR -> b^2 Ljmr
alpha = (Ljmr/Lq)^(1/4); 
Q = alpha*Q; 
J = J/alpha; 
R = R/alpha;  

disp('Display of iteration number and error ||A-(J-R)*Q||_F: ');  

i = 1; 
while i <= maxiter && cputime-cput <= timemax
    [Q] = porthamQfg(A,(J-R),Q); 
    [J,R] = porthamJRfg(A,J,R,Q); 

    e(i+1) = norm( A - (J-R)*Q , 'fro' ); 
    t(i+1) = cputime - cput;  
    
    fprintf('%2.0f:%2.3f - ',i,e(i+1)); 
    if mod(i,10) == 0
        fprintf('\n');       
    end
    i = i+1; 
end
fprintf('\n'); 