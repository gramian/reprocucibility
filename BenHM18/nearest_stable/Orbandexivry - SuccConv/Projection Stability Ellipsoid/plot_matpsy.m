%This script can be executed be the user while in debug mode in
%find_lambda.m 
%It produces a plot of the function psy for certain values of lambda.
%Find_lambda tries to find the point where ypsy is equal to zero. 

slab = lambda;
clear ypsy

lambda = [slab-0.9:0.001:slab];
nla = length(lambda);

for a = 1:nla
    ypsy(a) = EvalPsy(Xbar,Dp,Dq,lambda(a),n);
end

plot(lambda,ypsy)

lambda = slab