% Find the nearest stable matrix in the form (J-R)*Q to a given matrix A
% 
% ****** Description ******
% This code tries to solve the following optimization problem
% 
% min_{J,R,Q} || A - (J-R) * Q ||_F 
% 
% such that J = -J', R and Q PSD. 
%
% It uses a *fast gradient method* from smooth convex optimization,
% with a safety procedure. 
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

function [J,R,Q,e,t] = stableLinearFGM(A,maxiter,timemax,J,R,Q); 

cput = cputime; 
n = size(A,1); 
if nargin < 2
    maxiter = 1e6; 
end
if nargin < 3
    timemax = 10; 
end
if nargin < 4 || isempty(Q)
    Q = eye(n); 
    J = (A-A')/2; 
    R = projectPSD(J-A); 
    ejr =  norm( A - (J-R)*Q , 'fro'); 
end
% Scaling (J,R,Q) for better numerical stability
% -- For Dj and Dr
[uq,Lq,vq] = svds(Q,1); 
Lq = Lq^2; %=max(eig(Q*Q'))
% -- For Dq
JmR = J-R; 
JmRtJmR = JmR'*JmR; 
[ujmr,Ljmr,vjmr] = svds(JmRtJmR,1); 
% -- scaling such that Lq = Ljmr: a*Q -> a^2Lq et b b*JmR -> b^2 Ljmr
scal = (Ljmr/Lq)^(1/4); 
Q = scal*Q; 
J = J/scal; 
R = R/scal;  
L = scal^2*Lq; 

e(1) =  norm( A - (J-R)*Q , 'fro'); 
t(1) = cputime-cput;  
step = 1/L; 

i = 1; 
alpha0 = 0.1; % Parameter, can be tuned. 
alpha(1) = alpha0;

Yj = J; 
Yr = R; 
Yq = Q; 
restarti = 1; 

lsparam = 1.5; 
inneritermax = 10;  % lsparam^inneritermax approx 1e-6 
% Usually, when it is necessary to reduce the original step by a value
% larger than 1.5^10, it means there is no descent direction and it is
% better to reinitialize the FGM 

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

while i <= maxiter && cputime-cput <= timemax
    alpha(i+1) = ( sqrt(alpha(i)^4 + 4*alpha(i)^2 ) - alpha(i)^2) / (2);
    beta(i) = alpha(i)*(1-alpha(i))/(alpha(i)^2+alpha(i+1));
    
    D = (Yj-Yr)*Yq - A;
    gJ = D*Yq';
    gR = -gJ;
    gQ = (Yj-Yr)'*D;
    
    e(i+1) = +Inf;
    inneriter = 1;
    while e(i+1) > e(i) && inneriter <= inneritermax
        Jn = Yj - gJ*step;
        Jn = (Jn-Jn')/2;
        Rn = projectPSD(Yr - gR*step);
        Qn = projectPSD(Yq - gQ*step);
        
        e(i+1) =  norm( A - (Jn-Rn)*Qn , 'fro');
        t(i+1) = cputime-cput;
        
        step = step/lsparam;
        inneriter = inneriter+1;
    end
    
    if inneriter == inneritermax+1
        if restarti == 1
            % Restart if not a descent direction
            restarti = 0;
            alpha(i+1) = alpha0;
            Yj = J;
            Yr = R;
            Yq = Q;
            e(i+1) = e(i);
            % reduce lsparam
            %lsparam = max(1.01, (lsparam+1)/2);
            step = 1/L;
            disp('Descent could not be achieved: restart.')
        elseif restarti == 0 % no previous restart and no descent direction => converged to a stationary point
            disp('The algorithm has converged.')
            return;
        end
    else
        restarti = 1;
        % Conjugate
        Yj = Jn + beta(i)*(Jn-J);
        Yr = Rn + beta(i)*(Rn-R);
        Yq = Qn + beta(i)*(Qn-Q);
        
        % Scaling (J,R,Q) for better numerical stability and recomputing
        JmR = Jn-Rn;
        % Power method
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
        
        % -- scaling such that Lq = Ljmr: a*Q -> a^2Lq et b b*JmR -> b^2 Ljmr
        scal = (Ljmr/Lq)^(1/4);
        
        Qn = scal*Qn;
        Jn = Jn/scal;
        Rn = Rn/scal;
        
        L = (scal^2*Lq + Ljmr/scal^2)/2;
        
        Yj = Yj/scal;
        Yr = Yr/scal;
        Yq = Yq*scal;
        
        step = 1/L;
        
        J = Jn;
        R = Rn;
        Q = Qn;
    end
    
    if mod(i,ndisplay) == 0
        fprintf('%2.0f:%2.3f - ',i,e(i+1));
    end
    if mod(i,ndisplay*10) == 0
        fprintf('\n');
    end
    
    i = i+1;
end
fprintf('\n');