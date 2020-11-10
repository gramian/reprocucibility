% Fast gradient method to solve 
% solve min_{J,R} || A - (J-R) Q || such that J=-J', R PSD 
% 
% ****** Input ******
% A, Q   : two square matrices 
% (J,R)  : initialialization for variables (J,R) 
% maxiter: max number of iterations (default = 1e4)
% ****** Output ******
% (J,R)   : PSD matrix R, J=-J' such that || A - (J-R) Q || is minimized 
% e       : evolution of the error. 

function [J,R,e] = porthamJRfg(A,J,R,Q,maxiter,delta) 
 
n = size(A,1); 
if nargin <= 4
    maxiter = 1000;
end
if nargin <= 5
    delta = 1e-2; 
end
% Hessians, Lipschitz constants and scaling 
QQt = Q*Q'; 
[u,L,v] = svds(QQt,1); 

alpha0 = 0.1; % Parameter, can be tuned. 
alpha(1) = alpha0;
Yj = J; 
Yr = R; 

% Precomputation
AQt = A*Q'; 
i = 1; delta = 1e-3; 
eps = 1; eps0 = 0; 

% Fast gradient 
while i <= maxiter && eps >= delta*eps0
    % Previous iterates
    Jp = J; Rp = R; 
    % FGM Coefficients  
    alpha(i+1) = ( sqrt(alpha(i)^4 + 4*alpha(i)^2 ) - alpha(i)^2) / (2); 
    beta(i) = alpha(i)*(1-alpha(i))/(alpha(i)^2+alpha(i+1)); 
     
    grad = AQt - (J-R)*QQt; 
    
    R = projectPSD( Yr - grad / L  );
    
    J =  Yj + grad / L  ;
    J = (J-J')/2; 
    
    % `Optimal' linear combination of iterates
    Yr = R + beta(i)*(R-Rp); 
    Yj = J + beta(i)*(J-Jp); 
    
    % Error
    if nargout >= 3 
        e(i) = norm( A - (J-R)*Q , 'fro' ); 
    end
 
    % Evolution of iterates
     if i == 1
         eps0 = norm(Jp-J,'fro')/(norm(J,'fro')+1e-6) + norm(Rp-R,'fro')/(norm(R,'fro')+1e-6); 
     end
     eps = norm(Jp-J,'fro')/(norm(J,'fro')+1e-6) + norm(Rp-R,'fro')/(norm(R,'fro')+1e-6); 
     
     i = i + 1; 
end  