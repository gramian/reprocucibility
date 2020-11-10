function [f,g] = diststable(x, pars)
% function and gradient of ||X-A||_F + rho*max(0, alpha(X))
% where alpha denotes spectral abscissa and rho is a penalty parameter
% gradient is defined almost everywhere
A = pars.A;
if isreal(A)
    complx = 0;
else % in this case optimization is in space of 2n^2 real variables
    complx = 1;
end
penalty = pars.penalty;
if complx == 0
    n = length(A);
    N = n^2;
    X = reshape(x,n,n);
else
    n = length(A);
    N = 2*n^2;
    X = reshape(x(1:N/2),n,n) + 1i*reshape(x(N/2+1:N),n,n);
end
% f = 0.5*norm(X-A,'fro')^2; % using squared norm is not so effective
% g = reshape(X-A,n*n,1);
f = norm(X-A,'fro');
G = (X-A)/f;
[V,Lambda] = eig(X);
lambda = diag(Lambda);
% do not worry about breaking ties -- see papers for justification
[absc,indx] = max(real(lambda));
v = V(:,indx); % corresp right evector (column)
W = inv(V); % matrix of left eigenvectors
% pert = 1e-16;
% while isnaninf(W)
%     % perturb V until it is numerically nonsingular
%     V = V + pert*randn(n,n);
%     W = inv(V);
%     pert = pert*10
% end 
w = W(indx,:); % corresp left evector
if  absc >= 0
    f = f + penalty*absc;
    % gradient of abscissa: see Burke and Overton, 2001 for details
    G = G + penalty*w'*v'; % outer product of conjugate transposes of left and right eigenvector
end
if complx == 0
    g = real(reshape(G,n^2,1)); % discard imaginary rounding error
else
    % some thought is needed to justify this.  The penalty function has two
    % terms, both real functions of complex variables: the distance term
    % and the abscissa term. Let phi denote either term and let Grad
    % denote its gradient in complex matrix space.  Then we have
    %    phi(X + dX) - phi(X) \approx <Grad,dX> where <,> is the real trace
    % inner product on complex matrix space, so
    %    phi(X + dX) - phi(X) \approx Re trace Grad^* dX
    %  = (vec Re Grad)^T (vec Re dX) + (vec -i*Im Grad)^T (vec i*Im dX)
    %  = (vec Re Grad)^T (vec Re dX) + (vec Im Grad)^T (vec Im dX)
    % so setting g to the real and imaginary parts of G is the right thing.
    g = [reshape(real(G),n^2,1); reshape(imag(G),n^2,1)];
end