function pars = pars_diststable(n,penalty,complx)
% define "pars" for distance to stable matrix penalty function
if nargin < 3 
    complx = 0;
end
pars.fgname = 'diststable';
if complx == 0
    pars.nvar = n*n;
    pars.A = randn(n);
else
    pars.nvar = 2*n*n;
    pars.A = randn(n) + i*randn(n);
end
if nargin < 2
    pars.penalty = 1;
else
    pars.penalty = penalty;
end