function [psy] = PolyEvalPsy(ax,F,G,lambda,ellips_radius,n)
%Evaluate the function psy (formula 21 on page 6 of the notes)
%This function is used in PolyFindLambda

%History : 28/09/2010: Creation
%Last modified : 05/01/2011 : Cleaning of the file

%numerator
ax2 = ax.^2;
num =F*ax2;

%denominator
den = (diag(G)+lambda*diag(F)).^2;

%numerator over denominator
Mdiag = num./den;

%we perform a frobenius scalar product
psy = real(sum(Mdiag))-ellips_radius;

