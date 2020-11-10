% Sovles 
% 
% P solves: 
% P PSD and PA + A*P <= 0 


function P = FindstableP(X,param) 

n = size(X,1); 

if nargin <= 1
    param = 1e-16;
end

cvx_begin quiet sdp 
    variable P(n,n);
    variable lambda(1) 

    maximize( lambda   );

    subject to
    
        (P - param*eye(n) ) == hermitian_semidefinite( n ); 
        
        (-X*P - P*X' - lambda*eye(n)) == hermitian_semidefinite( n ); 
        
        norm(P,'fro') <= sqrt(n); 
        
        %lambda >= 0; 
        
cvx_end

if isnan(P(1,1))
    disp('The matrix may not be stable, the SDP could not be solved'); 
end

egpxxp = eig((-X*P - P*X')); 
if min( egpxxp ) < 0
    fprintf('An eigenvalue of -XP - PX* is negative %2.2d', min(egpxxp));
end