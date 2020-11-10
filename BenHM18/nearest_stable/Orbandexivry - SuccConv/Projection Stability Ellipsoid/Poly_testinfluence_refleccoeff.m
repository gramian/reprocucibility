%This routine is aimed at showing the behaviour of the roots of the
%starting point with respect to the way we modify the "instable" reflection
%coefficients (those that are greater than one). 
%In order to achieve our goal, 'PolySPReflection' need to be modified befor
%each run. 
%Three runs were initially achieved:
%       - 'k linear' : The instable reflection coefficients were replaced
%                      by a simple coefficient. k = coeff;
%       - 'k proportional' : The inst. refl. coef. were replaced by a
%                            reflection coef. propotional to themselves,
%                            aka : k = coeff*k;
%       -'k  inv. prop.' : The inst. refl. coef. were replaced by a refl.
%                          coef. inversely proportional to them, aka : 
%                           k = coeff*k/abs(k)^2 
%                          Note that k/abs(k)^2 is the mirrored refl. coef.
%                          with respect to the unit circle.



%Original polynomial
PolyStableData();

%defining the coeff.
k = linspace(-1,1,500);
len = length(k);

n = size(A,2);
a = poly(A);

%roots of the stable starting point 
root = zeros(len+1,n);
root(1,:) = roots(poly(A));

%distance to the original instable polynomial
objval = zeros(len,1);

%finding the starting point for each modification of the coeff
for i = 1:len
    b0 = PolySPReflection(A,e_lon,isdiscrete,k(i));
    b0 = [1 fliplr(b0.')];
    root(i+1,:) = roots(b0);
    objval(i) = norm(a-b0,'fro');
end

%plotting the path of the eigenvalues of the coeff
close all;
figure;
plot_iterates(root,isdiscrete,e_lon,roc,savf,path)
figure;
plot(k,objval)
ylabel('Distance to a(z)')
    
    
