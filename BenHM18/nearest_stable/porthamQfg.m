% Fast gradient method to solve 
% solve min_{Q} || A - U Q || such that Q PSD 
% 
% ****** Input ******
% A, U    : square matrices 
% Q       : initialialization for variable Q 
% maxiter : max number of iterations (default = 1e4)
% ****** Output ******
% Q       : PSD matrix such that || A - U Q || is minimized 
% e       : evolution of the error. 

function [Q,e] = porthamQfg(A,U,Q0,maxiter,delta)

n = size(A,1); 

if nargin <= 3
    maxiter = 1000;
end
if nargin <= 4
    delta = 1e-2;
end
if nargin <= 2 || isempty(Q0) 
    Q0 = eye(n); 
end

nA = norm(A,'fro')^2; 
Q = Q0; 

% Hessian and Lipschitz constant 
UtU = U'*U; 
[u,L,v] = svds(UtU,1); 
% Linear term 
UtA = U'*A; 

alpha0 = 0.1; % Parameter, can be tuned. 
alpha(1) = alpha0;
Y = Q; 
i = 1; 

% Stop if ||V^{k}-V^{k+1}||_F <= delta * ||V^{0}-V^{1}||_F
eps0 = 0; eps = 1;  

ebest = norm(A - U*Q , 'fro'); Qbest = Q;  

while i <= maxiter && eps >= delta*eps0
    % Previous iterate
    Qp = Q; 
    
    % FGM Coefficients  
    alpha(i+1) = ( sqrt(alpha(i)^4 + 4*alpha(i)^2 ) - alpha(i)^2) / (2); 
    beta(i) = alpha(i)*(1-alpha(i))/(alpha(i)^2+alpha(i+1)); 
    
    % Projected gradient step from Y
    Q = projectPSD( Y - (UtU*Y-UtA) / L );
    
    % `Optimal' linear combination of iterates
    Y = Q + beta(i)*(Q-Qp); 
    
    % Error 
    if nargout >= 2
        e(i) = norm(A - U*Q , 'fro');  
    end
    
    if i == 1
        eps0 = norm(Q-Qp,'fro'); 
    end
    eps = norm(Q-Qp,'fro'); 
    
%     % FG with restart: it may work better  
%     if i >= 2 && e(i) > e(i-1)
%        Y = Q; 
%        alpha(i+1) = alpha0;
%     end
    
    i = i + 1; 
end 