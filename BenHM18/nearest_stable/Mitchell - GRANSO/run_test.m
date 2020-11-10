% TEST 
clear all; clc; 

addpath('./granso');

n = 20; 
A = grcar(n,3); 

save A A; 

% Initialization 
N = n^2;
J = (A-A')/2; 
R = projectPSD(J-A); 
x0 = reshape(J-R,N,1);

compute_all_fn = makeFunctionsCombined();

% Options of GRANSO 
options.x0 = x0; 
options.maxit = 500; 
options.maxclocktime = 300; 
options.quadprog_info_msg = false; 
options.opt_tol = 1e-12; 

tic; 
soln2 = granso(N,compute_all_fn,options); 
toc; 
X2 = reshape(soln2.final.x,n,n); 
fprintf('Final error GRANSO = %2.2f, feasibility: %d.\n', norm(A-X2,'fro'), max(real(eig(X2)))); 