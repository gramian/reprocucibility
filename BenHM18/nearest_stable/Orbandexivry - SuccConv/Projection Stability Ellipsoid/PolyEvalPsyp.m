function [psyp] = PolyEvalPsyp(ax,F,G,lambda,n)
%Evaluate the derivative of the function psy (formula 21 on page 6 of the
%notes)
%This function is used in PolyFindLambda

%History : 28/09/2010: Creation
%          05/01/2011: Cleaning of the file
%          05/06/2012: Introducing G
%Last modified : 05/01/2011 : Cleaning of the file

%numerator
ax2 = ax.^2;
num =F.^2*ax2;

%denominator
den = (diag(G)+lambda*diag(F)).^3;

%numerator over denominator
Mdiag = num./den;

%we perform a frobenius scalar product
psyp = -2*real(sum(Mdiag));
