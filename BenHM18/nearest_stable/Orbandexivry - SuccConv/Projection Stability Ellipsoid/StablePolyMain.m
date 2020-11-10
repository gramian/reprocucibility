function [bi,s]= StablePolyMain(a,s,isdiscrete,compute,verbose,e_lon,roc,savf,path)
%[bi]= StablePolyMain(a,S,isdiscrete,compute,verbose,e_lon,roc,savf,path)
%StablePolyMain find the closest stable polynomial to 'a' in the sense of 
%the frobenius norm without explicitly computing its roots. The algorithm
%transforms the polynomial into a companion matrix which is the sum of a
%shift matrix and a rank one matrix. The eigenvalues of the companion
%matrix are identical to the roots of the polynomial it represents.
%
%   [bi] = StablePolyMain(a,S,isdiscrete,compute,verbose) returns a
%          polynomial 'bi' of the same degree than a but with stable roots.
%          The stable roots of 'a' can be modified. The polynomial 'bi' is
%          computed as the local minimum of an optimization problem.
%
%          The algorithm stops when the change between two successive
%          iterations gets smaller than 1e-4(default). If this criterion
%          cannot be satisfied then the algorithm exits with an error.
%
%   Input: a : a monic polynomial of degree n, or a companion matrix of
%          size n by n representing a polynomial.The polynomial a is
%          obtained using the 'poly()' command.
%          isdiscrete : 1 if 'a' describes a discrete time system
%                       0 if 'a' describes a continuous time system
%          compute    : if 1 the algorithm will compute and keep track of
%          the eigenvalues of the system in the file stock.mat. In that
%          file the first line gives the original eigenvalues, the second
%          line gives the eigenvalues of the starting point matrix S , the
%          following lines gives the evolution of the eigenvalues. 
%          A graph showing the computed values is produced as well. (Blue
%          circles : starting point, red circles : initial eigenvalues,
%          darkness gives the evolution of the iterations). 
%          verbose : if 1 print information during the computation
%          S : a starting point, in companion form whose dimension is equal
%          to the degree of the polynomial. Its eigenvalues must be
%          stable.
%          e_lon : the "distance" (+ or -) from the stability boundary at
%          which we want the 'stabilized' eigenvalues to lie inside the
%          stability region. 
%       if compute == 1
%          roc  : if 1 gives a semilog plot of the convergence of each eigenvalue
%          savf : if 'None' will not save the graphs, otherwise will save
%          them after formatting in the the directory given by path. Other
%          values for savf are: 'Article', 'Poster', 'Personnalized' or
%          'Default'. Values linked to those keywords can be edited in
%          Stable_param.m.
%          path : the path to an existing directory.
%
%
%  Output: bi : a stable polynomial whose roots are the closest to a's.
%
%History : 05/01/2011 : Cleaning of the unused code from the file, a backup
%                       (in the bup file dated 31/12/2010) exists which
%                       features a history of the resolution method
%                       'Find_PH'.
%          ../03/2012 : A test is created to see if random choices of P
%                       and/or Q brings a comparable performance. For this
%                       the biggest change is the addition of a structure
%                       's' to keep track of the generated P and Q. Some
%                       lines must be actived/deactived for using this
%                       option.
%          24/08/2012 : A barrier term is added after a full backup of the
%                       various files. Everuthing is copied from the matrix
%                       case.
%Last modified : 24/08/2012

%% preprocessing


%adding the Utils directory to matlab path.
d = pwd;
smotd = length(d);
cd ..
str =pwd;
smots = length(str);
d = d(smots+2:smotd);
cd(d);
str = strcat(str,'\Utils');
addpath(genpath(str));


%preprocessing
if nargin < 2
    S = [];
else
    [nls ncs] = size(s);
    
    if nls < ncs
        %if polynomials are given as a line, they will be treated as if they
        %were obtained using the poly() command
        s = s.';
        s= flipud(s);%We work with column polynomials whose last element is 1
        temp = nls;
        nls = ncs;
        ncs = temp;
    end
    if nls == ncs
        S = s;
    elseif ncs == 1
        %the last index must be one
        firstelem = s(nls);
        if firstelem ~=1
            error('The polynomial ''s'' is not monic, please consider using the poly command')
        end
        
        s(nls) = [];
        nls = nls-1;
        S = diag(ones(nls-1,1),-1);
        S(:,nls) = -s;
    else
        error('-> Input "s" is not a vector or a square matrix representing a polynomial');
    end
end

if nargin < 3
    isdiscrete = 0;
end
if nargin < 4
    compute = 1;
end
if nargin < 5
    verbose = 0;
end
if nargin < 6
    e_lon = 0;
end
if nargin < 9
    roc  = 0;%rate of convergence
    savf = 'none';%type of save
    path = '';
end
namefig = '';
if nargin == 0
    try 
        PolyStableData()
    catch
        error('If called without argument, StablePolyMain calls a data file named "PolyStableData.m", this file is nonexistent')
    end      
end
if (nargin > 6 && nargin < 9 ) || nargin > 9
    error('-> Wrong number of input arguments')
end


if strcmpi(savf,'none')~= 1
    if ~isdir(path)
        error('The given directory does not exist')
    end
end

[nl nc]= size(a);

if nl < nc
    %if polynomials are given as a line, they will be treated as if they
    %were obtained using the poly() command
    a = a.';
    a= flipud(a);%We work with column polynomials whose last element is 1
    nl = nc;
    nc = 1;
end

if nl == nc
    A = a;
    a = -A(:,nc);
elseif nc == 1
    %the last index must be one
    firstelem = a(nl);
    if firstelem ~=1
        error('The given polynomial is not monic, please consider using the poly command')
    end
        
    %a = a/firstelem;
    %the last index is not used in a companion matrix.
    a(nl) = [];
    nl = nl-1;
else
    error('-> Input "a" is not a vector or a square matrix representing a polynomial');
end
n = nl;

%checking if the problem is real or not, pbisreal is used for the
%computation of the hessian of the barrier function.
if isreal(a)%&& isreal(bi)
    pbisreal = 1;
else
    pbisreal =0;
end


%% processing
tic

%parameters
NbItereMax = 1e6;
ellips_radius = 0.99;
beta = 1+1/sqrt(n);%%WARNING : En fonction des problèmes, les valeurs
                          %%de beta 1.01 et 1+1/n apportent des solutions
 %1.01:benchmark          %%différentes. En aucun cas l'une ou l'autre de
 %1+1/n :le reste         %%ces possibilités n'amènent à une solution
                          %%automatiquement optimale. Dans le cas convex,
                          %%1+1/sqrt(n) is optimal pour matrices => 1+1/n
                          %%pour polynomes
tol = n*1e-8;%don't put it smaller than 1e-10 otherwise it might be under machine precision in its search of lambda
rel_mu_min = n*1e-10;%should not be higher than tol
deltatolmu = n*1e1;
if tol < 1e-10
    tol = 1e-10;
    warning('Cannot set to small relative accuracy :relative accuracy set to 1e-10'); %#ok<WNTAG>
end

%initialization                          
itere = 0;                          
la = 1;

%constructing the companion matrix
Z = diag(ones(n-1,1),-1);
en = [zeros(n-1,1);1];
A = Z - a*en';

%saving the iterations
if compute
    stock = [];
    Bmat = [];
    [stock, Bmat] = stockIt(stock,Bmat,A);
    VecSqrtObjval = [];
end

%This is one input arguement for MatFindPfromQrandom. When the line is put
%in commment this arguement changes at each iterations so that Q is
%randomized. When the line is activated Q is the same for all the
%iterations.
% Lq = [];
% Lq = 1e-1*rand(n);


%taking a starting point b0 to initialize the algorithm
bi = PolyFindStartingPoint(A,e_lon,isdiscrete,S);
Bi = Z - bi*en';
if compute
    [stock, Bmat] = stockIt(stock,Bmat,Bi);
    VecSqrtObjval = [VecSqrtObjval; norm(A-Bi,'fro')];
end


%Central path parameter
lmin = abs(max(real(roots([1 fliplr(bi.')]))));
mu  =min(1/(sqrt(n*lmin^2)),sqrt(n*lmin^2));


%entering the loop that will iterate to find the best B
Y = Bi;
logdetY = -log(real(det(Y)));
sqrtObjval = norm(A-Bi,'fro');
DiffNormObjval = 0.5*sqrtObjval^2;%initial step between A and the starting point
DiffObjFunction = DiffNormObjval; %negating the logdet contribution ensures that it goes at least once in the while loop)
mumin = rel_mu_min*DiffNormObjval/logdetY;
if mu < mumin
    mumin = mu;
end
hasconverged = 0;
try
    while (~hasconverged&& (itere < NbItereMax))

        %counting iteration
        itere = itere + 1;

    %---%Find best P and Q
        %same version as the matrix case
        [Dp,Dq,T,invT] = MatFindPQ(Z-bi*en',n,e_lon,isdiscrete,verbose);
        %The next line should be actived if Q has to be random
%         [Dp,Dq,T,invT] = MatFindPfromQrandom(Bi,Lq,n,e_lon,isdiscrete,verbose);
        P = T*Dp*T';
        iQ = invT'*diag(diag(Dq.^(-1)))*invT;
    
    %---Computing the hessian of the barrier function
        [F,w] = PolyFindF(a,bi,P,iQ,mu,pbisreal,isdiscrete);
        if pbisreal
            sizeF = n;
        else
            sizeF = 2*n;
        end
        
    %---It is possible to compute the eigenvalues of F so
        %that the problem is easier.
        [U,Df] = eig(F);
        Dg = eye(sizeF)+mu*Df;

        %independent term
        ax = U'*w;

    %---%Find the best lambda. The calculus are not depedent on the kind of
    %system (discrete or continuous)
        [la] = PolyFindLambda(ax,Df,Dg,la,sizeF,tol/1000,ellips_radius,verbose);

    %---%Solving the main problem
        uh = ax./(diag(Dg)+la*diag(Df));
        h = U*uh;
        
        if ~pbisreal
            h = complex(h(1:n),h(n+1:2*n));
        end

    %---%Updating Bi
        bi = bi+h;%vector of coefficient of the Companion matrix
        
    %---%Computing the update of the stopping criterion
        sqrtObjvalprev = sqrtObjval;
        sqrtObjval = norm(a-bi,'fro');
        DiffNormObjval = 0.5*(sqrtObjvalprev^2-sqrtObjval^2);
        logdetYprev = logdetY;
        ZP = [zeros(1,n);P(1:n-1,:)];
        ZPPZ= ZP+ZP';
        pn= P(:,n);
        bip = bi*pn.';
        if ~isdiscrete
            Y = ZPPZ-(bip+bip');
            Y = -(Y+Y');
        else
            [U,D] = eig(P);
            d = sqrt(diag(D));
            Pd = U*diag(d)*U';
            Y = [Pd 2*Pd*(Z-bi*en'); zeros(n) Pd];
            Y = (Y+Y')/2;
        end
        logdetY = -log(real(det(Y)));    
        DiffObjFunction = abs(DiffNormObjval + mu*(logdetYprev-logdetY));
        if isnan(DiffObjFunction) 
            error('Bad conditionning')
        end
         %--Update of mu
        if DiffObjFunction < tol*deltatolmu
            if ((DiffObjFunction < tol) && (mu <= mumin))
                hasconverged = 1;
            end
            mu = max(mu/beta,mumin);
        end

    %---%saving the evolution of the solution
        if compute
            [stock,Bmat] = stockIt(stock,Bmat,bi);
            VecSqrtObjval = [VecSqrtObjval;sqrtObjval]; %#ok<AGROW>
        end

    %---%print information if necessary
        if verbose
            fprintf('\n-> Iteration %4i : Distance to the original matrix : %8.6e\n',itere,sqrtObjval)
            fprintf('                     Improvement of the obj. val.   : %3.2e\n',DiffObjFunction)
            fprintf('                     Norm of the direction h        : %3.2e\n',norm(h,'fro'))
            fprintf('                     Value of the barrier param.    : %3.2e\n',mu)
            %correction_bi = [0 fliplr(h.')]%#ok<NASGU,NOPRT>
            %coeff_bi = [1 fliplr(bi.')] %#ok<NASGU,NOPRT>
        end
    end
    
    %returning a polynomial
    bi = [1 fliplr(bi.')]; %compatible with Matlab definition of polynomial
    bi = bi*firstelem;
    %opening figure for the evolution of the solution
    if compute
        close all
        plot_iterates(stock,isdiscrete,e_lon,roc,savf,path,namefig);
        VecSqrtObjval = [min(VecSqrtObjval);VecSqrtObjval];
        save stock.mat stock Bmat VecSqrtObjval
    end

    %displaying time
    time = toc;
    if verbose
        time
    end

catch
    toc
    %opening figure for the evolution of the solution
    err = lasterror;
    mess = err.message
    id = err.identifier
    stack = err.stack(1,1)
    err2throw = struct('identifier',id,'stack',stack,'message',mess);
    if compute
        if (strcmp(id,'control:InputArguments')) && (itere>1)
            stock(end,:) = [];
            Bmat(:,end) = [];
            VecSqrtObjval(end) = [];
        end
        close all
        plot_iterates(stock,isdiscrete,e_lon,roc,savf,path,namefig)
        VecSqrtObjval = [min(VecSqrtObjval);VecSqrtObjval];
        save stock.mat stock Bmat VecSqrtObjval
        display('Evolution of the Bi saved in stock.mat')
    end
    fprintf('Failure at iteration %i\n returns last stable result',itere)
    bi = bi - h;
    bi = [1 fliplr(bi.')]; %compatible with Matlab definition of polynomial
    bi = bi*firstelem;
    return
    %error(err2throw)
end

