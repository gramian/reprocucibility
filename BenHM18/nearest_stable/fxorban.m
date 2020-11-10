% Matrix from the paper 
% F.-X. Orbandexivry, Yu. Nesterov, and P. Van Dooren, Nearest stable 
% system using successive convex approximations, Automatica, 49 (2013), 
% pp. 1195-1203. 

function A = fxorban(n,delta) 

A = diag(ones(n-1,1),-1); 
A(1,n) = delta; 