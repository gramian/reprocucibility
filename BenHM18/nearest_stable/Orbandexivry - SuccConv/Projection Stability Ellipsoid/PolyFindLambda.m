function [lambda] = PolyFindLambda(ax,F,G,lambda,n,Tol,ellips_radius,verbose)
%[lambda] = PolyFindLambda(ax,F,lambda,n,Tol,verbose) computes the lagrange
%multiplier that solves a one dimensional equation found in PolyEvalPsy. It
%applies a newton scheme along with different switch to control
%convergence.
%
%Input : ax     : a parameter for EvalPsy (a-x)
%        F      : a parameter for EvalPsy
%        G      : a parameter for EvalPsy
%        lambda : the last known solution for lambda (1 if unknown)
%        n      : the degree of the polynomial a
%        Tol    : A precision to be achieved by the scheme
%        verbose: if 1, will print some information when needed.
%
%Output :lambda : the value of the lagrange multiplier solution of Psy.
%
%History : 28/09/2010: Creation
%          16/11/2010 : graphical check moved at the end of the method.
%          05/01/2011 : Making of the header
%          18/07/2011 : security on the value of Tol added
%          05/06/2012 : Remove maximal number of iteration, introducing G
%          29/12/2012 : Remove the minimum for the tol so put back the
%                       maximal number of iteration
%last modified : 18/07/2011

NbIterMax = 3000;


if (ellips_radius >=1)
    error('FindLambda:TooBigEllipsRadius','The radius of the ellipsoid is set greater or equal to 1 in Stable_param.m')
end
% %security on the tol (a tol at 1e-16 would created an infinite loop)
% if Tol < 1e-12
%     Tol = 1e-12;
% end

% if lambda < 5e-3%I know this is weird but starting form la_old/50 is good until la gets small.
%     lambda = lambda*1e1;%After the starting point for la, should be taken large enough.
%     %To be usefull those factor should be linked to the size of A.
% end
lambda_start = lambda;%conserving the first value of lambda

iter = 0;

%iterating on lambda
psy = 10;%the value given here aims at passing the while criterion
psyp = 1;
wquit = 0; %enables an early exit from the loop
while ((abs(-psy/psyp) > Tol)&& ~wquit )&&(iter < NbIterMax)

    %evaluating the function and its derivative
    psy = PolyEvalPsy(ax,F,G,lambda,ellips_radius,n);
    psyp= PolyEvalPsyp(ax,F,G,lambda,n);
    if abs(psyp) < eps
        if verbose
            warning('    The derivative of psy is zero : divergence of the newton method observed') %#ok<WNTAG>
        end
        %if abs(lambda)< eps %this case should only occur if the matrix A is already stable instead of being instable
            lambda =0;
        %end
        return
    end

    %actualizing the value of lambda
    lambda = lambda - psy/psyp;

    %exiting if lambda gets negative
    if lambda < 0
        lambda = lambda_start/10;%we start again the problem with a smaller lambda
        lambda_start = lambda;
        if lambda_start < eps % if the root of Psy is indeed negative
            wquit = 1;
        end
    end

    %actualizing the number of iteration
    iter = iter+1;
end

%graphical check of optimality
% Plot_polypsy %for debug purpose only

%warning generated if the maximum number of iteration is reached
if verbose && (iter >= NbIterMax)
    warning('-> Maximal number of iteration(%i) reached in FindLambda.m. At the last iteration the change in lambda was still of %5.3e',NbIterMax,-psy/psyp); %#ok<WNTAG>
end



