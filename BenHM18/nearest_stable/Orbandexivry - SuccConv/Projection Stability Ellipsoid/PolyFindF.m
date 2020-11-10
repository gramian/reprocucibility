function [F,y] = PolyFindF(a,bi,P,iQ,mu,pbisreal,isdiscrete)
%PolyFindF(Bi,P,iQ,isdiscrete,,pbisreal) computes the expression of the
%hessian of the barrier function based on its input. The expression of the
%hessian varies according to the type of the studied system, whether
%continuous or discrete.
%
% Input  : pbisreal : is 1 if the problem is real so that h has no
%                     imaginary part.
%
% Output : F : The hessian of the barrier function
%
%History : 15/11/2010
%          02/12/2010 : Correction of the sign of p(n), the good sign
%                             is a plus. (line 37).
%          18/07/2011 : change of F = 2*F in F = F+ F'
%          24/08/2012 : Insertion of a barrier function in the objective
%                       function so that convergence can be guaranteed. Now
%                       the independent term y of the equation F h = y is
%                       computed here as well
%last modified : 24/08/2012

n = length(bi);
Z = diag(ones(n-1,1),-1);

X = Z-[zeros(n,n-1) bi];
pn = P(:,n);
    
if ~isdiscrete
    iQp = iQ*pn;

    %t_1 = coefficient for h
    temp = real((pn'*iQp));
    t1 = temp*iQ;
    %t_2 = coefficient for conj(h)
    t2 = iQp*iQp.';
    
    %additional term independent of h
    addterm = 2*mu*iQp;
else
    Xp = X*pn;
    iQXp = iQ*Xp;
    
    %see notes on polynomial for formulae for t_1 and t_2
    %coefficient for h
    t1 = iQ*((Xp'*iQXp)+pn(n));
    %coefficient for conj(h)
    t2 = iQXp*iQXp.';
    
    %additional term independent of h
    addterm =  2*mu*iQXp;
end

t1re = real(t1);
t1im = imag(t1);
t2re = real(t2);
t2im = imag(t2);

%independent term in the equation F h = y
ytemp = a-bi+addterm;
if ~pbisreal
    %hessian of the barrier
    F = [t1re+t2re -t1im+t2im; ...%variable : [h_re;h_im];
         t1im+t2im t1re-t2re];
    y = [real(ytemp);imag(ytemp)];
else
    %if the problem is purely real
    F = t1re+t2re;
    y = ytemp;   
end

F = F+F';%(=2*F) the constant 2 is part of the hessian F, but at the same F must be perfectly symmetric