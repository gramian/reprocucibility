function [F] = PolyFindF(Bi,P,iQ,pbisreal,isdiscrete)
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
%Created       : 15/11/2010
%Modified      : 02/12/2010 : Correction of the sign of p(n), the good sign
%                             is a plus. (line 37).
%                18/07/2011 : change of F = 2*F in F = F+ F'
%Last modified : 18/07/2011

%P is square
n = size(P,1);
en = [zeros(n-1,1);1];

X = Bi;
p = P*en;
iQp = iQ*p;

if ~isdiscrete
    %t_1 = coefficient for h
    temp = real((p'*iQp));
    t1 = temp*iQ;
    %t_2 = coefficient for conj(h)
    t2 = iQp*iQp.';
else
    Xp = X*p;
    iQXp = iQ*Xp;
    
    %see notes on polynomial for formulae for t_1 and t_2
    %coefficient for h
    t1 = iQ*((Xp'*iQXp)+p(n));
    %coefficient for conj(h)
    t2 = iQXp*iQXp.';
end

t1re = real(t1);
t1im = imag(t1);
t2re = real(t2);
t2im = imag(t2);

if ~pbisreal
    %hessian of the barrier
    F = [t1re+t2re -t1im+t2im; ...%variable : [h_re;h_im];
         t1im+t2im t1re-t2re];
else
    %if the problem is purely real
    F = t1re+t2re;
end

F = F+F';%(=2*F) the constant 2 is part of the hessian F, but at the same F must be perfectly symmetric