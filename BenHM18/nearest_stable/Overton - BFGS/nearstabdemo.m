function [X,erec,trec] = nearstabdemo(A, maxiter, timemax, rho_init, stab_tol)
% Michael Overton, July 4, 2012 (complex version added Sept 25, 2012)
% call HANSO (Hybrid Algorithm for NonSmooth Optimization) to try to find
% the nearest stable matrix X to a given matrix A by minimizing the
% nonsmooth, non-Lipschitz penalty function
%                   ||X-A||_F + rho*alpha(X)  if X is unstable
%                   ||X-A||_F  otherwise
% where alpha is the spectral abscissa, max(real(eig(X))), and
% rho is a penalty parameter.  We start with a specified initial
% value for rho and repeatedly multiply by 10 if necessary.

addpath('.\hanso2_1'); 

if nargin < 1
    n = input('A will be generated randomly.  What size? ');
    complx = input('enter 0 for real, 1 for complex ');
    if complx == 0
        A = randn(n);
    else
        A = randn(n) + 1i*randn(n);
    end
else
    n = length(A);
    complx = ~isreal(A);
end
if nargin < 2
    maxiter = 1000; 
end
if nargin < 3
    timemax = 60; 
end
if nargin < 4
    rho_init = n^2; % initial value of penalty parameter (heuristic)
end
if nargin < 5 % tolerance for declaring X to be stable
              % cannot make too tight as multiple evalues v sensitive
    stab_tol = 1e-3;
end  
pars.A = A;
pars.fgname = 'diststable';
% use vec(A) as one of many starting points, but need to perturb it so
% objection function is differentiable

% if complx == 0
%     N = n^2;
%     x0 = A(:) +0.01*randn(N,1);
% else
%     N = 2*n^2;
%     x0 = [real(A(:)) + 0.01*randn(N/2,1); imag(A(:)) + 0.01*randn(N/2,1)];
% end
if complx == 0
    N = n^2;
    J = (A-A')/2; 
    R = projectPSD(J-A); 
    
    x0 = reshape(J-R,N,1); 
else
    N = 2*n^2; 
    J = (A-A')/2; 
    R = projectPSD(J-A); 
    X0 = (J-R); 
    x0 = [real(X0(:)) ; imag(X0(:)) ];
end


pars.nvar = N;
% augment this starting point with other random starting points
%x0 = [x0 randn(N,10)];

pars.penalty = rho_init;
if max(real(eig(A))) < 0
    fprintf('Nothing to do: A is stable\n')
    X = A;
    return
end
options.x0 = x0;
options.maxit = maxiter;
options.cpumax = timemax;
options.samprad = []; 
 if N > 1000 
     options.nvec = 50; 
 end

warning off % to avoid a lot of warnings about nearly singular systems
[x, f, ~, ~, ~, ~, ~, erec, trec] = hanso(pars,options);  % look for local minimizer
%[x, f] = hanso(pars,options); 
erec = erec{1}; 
trec = trec{1}; 

% while trec(end) < timemax
%     options.x0 = x;
%     options.maxit = maxiter-length(erec);
%     options.cpumax = timemax-trec(end);
%     [x, f, ~, ~, ~, ~, ~, erec2, trec2] = hanso(pars,options);
%     erec = [erec erec2{1}]; 
%     trec = [trec trec2{1}+trec(end)]; 
% end

if complx == 0
    X = reshape(x,n,n); % turn vector x into matrix X
else
    X = reshape(x(1:N/2),n,n) + 1i*reshape(x(N/2+1:N),n,n);
end
eigvalsX  = eig(X);
while max(real(eigvalsX)) >= stab_tol  % stability tolerance
    %if pars.penalty > 1e6
    fprintf('**************************************************************************'); 
    fprintf('**************************************************************************'); 
        fprintf('!!! It could not find a stable matrix: try better starting point !!!')
      fprintf('**************************************************************************'); 
      fprintf('**************************************************************************'); 
        %return
    %end
%     pars.penalty = 10*pars.penalty;
%     fprintf('\n best X found is not stable: increase penalty parameter to %g\n', pars.penalty)
%     pars.x0 = x; % use previous result as starting point
%     [x,f] = hanso(pars,options);
%     if complx == 0
%         X = reshape(x,n,n); % turn vector x into matrix X
%     else
%         X = reshape(x(1:N/2),n,n) + 1i*reshape(x(N/2+1:N),n,n);
%     end
%     eigvalsX = eig(X);
end 
%format long e
fprintf('\nDistance from A to Best X Found is %g\n', norm(X-A,'fro'))
fprintf('\nEigenvalues of Best X Found Are\n')
eigvalsX
[evecX,junk] = eig(X);
fprintf('Singular Values of its Eigenvector Matrix Are\n')
singvalsV = svd(evecX)
fprintf('Small singular values indicate that exactly optimal X may have multiple eigenvalues.\n')
fprintf('In principle, local optimality conditions for X can be checked by approach outlined in\n')
fprintf('Theorem 10.2, Variational Analysis of Non-Lipschitz Spectral Functions, Burke&Overton, 2001.\n')
fprintf('However, this requires computing the Jordan form of X which is nontrivial.\n')
evalA = eig(pars.A);
% % plot eigenvalues of original A and proposed optimal X
% figure(1)
% clf
% plot(real(evalA),imag(evalA),'bx');
% hold on
% plot(real(eigvalsX),imag(eigvalsX),'ro');
% legend('eigenvalues of original A','eigenvalue of proposed X')
% title('Eigenvalues of Original A and Candidate Nearest Stable Matrix X')
