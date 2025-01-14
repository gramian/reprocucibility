%This script can be executed be the user while in debug mode in
%PolyFind_lambda.m 
%It produces a plot of the function psy for certain values of lambda.
%Find_lambda tries to find the point where ypsy is equal to zero. 

slab = lambda;
clear ypsy
interv = slab;
if slab >10
    interv = 10;
end
    
% lambda = [slab-0.9*interv:0.001:slab+0.9*interv];
lambda = [0.001:0.1:slab+10];
% lambda = linspace(1e-4,slab+10,1000);
% lambda = [0.999*slab:1e-4*slab:1.001*slab];
nla = length(lambda);
ypsy = zeros(1,nla);

for a = 1:nla
    ypsy(a) = PolyEvalPsy(ax,F,G,lambda(a),ellips_radius,n);
end

plot(lambda,ypsy)
hold on;
xrange = get(gca,'xlim');
idline = line(xrange,[0 0]);
set(idline,'color','r')
hold off;

lambda = slab;