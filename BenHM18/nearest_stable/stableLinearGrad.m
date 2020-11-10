% Find the nearest stable matrix in the form (J-R)*Q to a given matrix A
% 
% ****** Description ******
% This code tries to solve the following optimization problem
% 
% min_{J,R,Q} || A - (J-R) * Q ||_F 
% 
% such that J = -J', R and Q PSD. 
%
% It uses a *projected gradient method*. 
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

function [J,R,Q,e,t] = stableLinearGrad(A,maxiter,timemax,J,R,Q); 

cput = cputime; 
% Initialization 
n = size(A,1); 
if nargin <= 1
    maxiter = 100; 
end
if nargin <= 2
    timemax = +Inf; 
end

if nargin <= 4 || isempty(Q)
    Q = eye(n); 
    J = (A-A')/2; 
    R = projectPSD(J-A); 
    ejr =  norm( A - (J-R)*Q , 'fro'); 
end

e(1) =  norm( A - (J-R)*Q , 'fro'); 
t(1) = cputime-cput;  

% Scaling (J,R,Q) for better numerical stability
[uq,Lq,vq] = svds(Q,1); 
Lq = Lq^2;  
JmR = J-R; 
JmRtJmR = JmR'*JmR; 
[ujmr,Ljmr,vjmr] = svds(JmRtJmR,1);  
alpha = (Ljmr/Lq)^(1/4); 
Q = alpha*Q; 
J = J/alpha; 
R = R/alpha;  
L = alpha^2*Lq; 

disp('Display of iteration number and error ||A-(J-R)*Q||_F: ');   
if n <= 50
    ndisplay = 1000;
elseif n <= 100
    ndisplay = 100;
elseif n <= 200
    ndisplay = 10;
else
    ndisplay = 1;
end

% parameter for the line search
lsparam = 2; 
step = 1/L; 

i = 1;  
while i <= maxiter && cputime-cput <= timemax
    D = (J-R)*Q - A; 
    gJ = D*Q'; 
    gR = -gJ; 
    gQ = (J-R)'*D; 
    
    e(i+1) = +Inf; 
    kin = 0; 
    while e(i+1) > e(i) && lsparam(1)^kin <= 1e16 
        Jn = J - gJ*step; 
        Jn = (Jn-Jn')/2; 
        Rn = projectPSD(R - gR*step); 
        Qn = projectPSD(Q - gQ*step); 
        
        e(i+1) =  norm( A - (Jn-Rn)*Qn , 'fro'); 
        t(i+1) = cputime-cput; 
        
        step = step/lsparam;
        kin = kin+1; 
    end
    if lsparam(1)^kin >= 1e16   
       disp('The algorithm has converged.') 
       return; 
    end
    JmR = Jn-Rn;
    % Scaling, using the Power method
    for powi = 1 : 4
        uq = Qn*uq; uq = Qn*uq; uq = uq/norm(uq);
        ujmr = JmR*ujmr; ujmr = JmR'*ujmr; ujmr = ujmr/norm(ujmr);
    end
    uq = Qn*uq;
    Lq = norm(uq)^2;
    uq = uq/norm(uq);
    ujmr = JmR*ujmr; ujmr = JmR'*ujmr;
    Ljmr = norm(ujmr);
    ujmr = ujmr/norm(ujmr);
    
    scal = (Ljmr/Lq)^(1/4);
    
    Qn = scal*Qn;
    Jn = Jn/scal;
    Rn = Rn/scal;
    
    L = (scal^2*Lq + Ljmr/scal^2)/2;
    
    step = 1/L; 
    
    J = Jn;
    R = Rn;
    Q = Qn;
    
    if mod(i,ndisplay) == 0
        fprintf('%2.0f:%2.3f - ',i,e(i+1)); 
    end
    if mod(i,ndisplay*10) == 0
        fprintf('\n');       
    end
    i = i+1; 
end
fprintf('\n');  