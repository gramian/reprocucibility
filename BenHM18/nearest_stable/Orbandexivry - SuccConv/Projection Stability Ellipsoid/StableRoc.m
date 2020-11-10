function  [convergence] = StableRoc()
%StableRoc compute the rate of convergence of the solution. The solution is
%considered as the last iterate obtained.
%Created: 23/09/2010 
%Lastmodified: 05/01/2010 : Comment
%              24/09/2010 : code

load stock.mat Bmat

[q,k] = size(Bmat);
n = sqrt(q);


A = Bmat(:,1);
Sol = Bmat(:,k);

A = reshape(A,n,n);
Sol = reshape(Sol,n,n);

bestnorm = norm(Sol-A,'fro');
convergence = zeros(k-2,1);

for i = 2 : k-1
    iterate = Bmat(:,i);
    itreshape = reshape(iterate,n,n);
    itnorm = norm(itreshape-A,'fro');
    convergence(i-1) = itnorm-bestnorm;
end


    


