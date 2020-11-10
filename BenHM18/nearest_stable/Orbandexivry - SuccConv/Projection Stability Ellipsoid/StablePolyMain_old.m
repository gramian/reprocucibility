function [bi,s]= StablePolyMain_old(a,S,isdiscrete,compute,verbose,e_lon,roc,savf,path)
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
%          size n by n representing a polynomial.
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
%
%Last modified : 05/01/2011

%% preprocessing
%for Nicolas M.
global s
s = struct('X',{},'P',{},'iQ',{});


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
end
if nargin < 3
    isdiscrete = 0;
end
if nargin < 4
    compute = 1;
end
if nargin < 5
    e_lon = 0;
end
if nargin < 8
    roc  = 0;%rate of convergence
    savf = 'none';%type of save
    path = '';
end
if nargin < 9
    verbose = 0;
end
if nargin == 0
    try 
        PolyStableData()
    catch
        error('If called without argument, StablePolyMain calls a data file named "PolyStableData.m", this file is nonexistent')
    end      
end
if (nargin > 5 && nargin < 8 ) || nargin > 9
    error('-> Wrong number of input arguments')
end


if strcmpi(savf,'none')~= 1
    if ~isdir(path)
        error('The given directory does not exist')
    end
end

nl = size(a,1);
nc = size(a,2);

if nl < nc
    a = a.';
    temp = nl;
    nl = nc;
    nc = temp;
end

if nl == nc
    A = a;
    a = -A(:,nc);
elseif nc == 1
    %the last index must be one
    if a(nl) ~=1
        a = a/a(nl);
    end
    %the last index is not used in a companion matrix.
    a(nl) = [];
    nl = nl-1;
else
    error('-> Input "a" is not a vector or a square matrix representing a polynomial');
end

tic

%  The problem is normed
n = nl;
Z = [zeros(1,n-1) 0;eye(n-1) zeros(n-1,1)];
en = [zeros(n-1,1);1];
A = Z - a*en';
itere = 0;
la = 1;
NbItereMax = 1e6;
tol = 1e-5;%don't put it smaller than 1e-10 otherwise it might be under machine precision in its search of lambda
if tol < 1e-10
    tol = 1e-10;
    warning('Cannot set to small relative accuracy :relative accuracy set to 1e-10');
end

%saving the iterations
if compute
    %for Nicolas M.
    s(1).X = A;
    
    stock = [];
    Bmat = [];
    [stock, Bmat] = stockIt(stock,Bmat,A);
    objval = [];
end

%This is one input arguement for MatFindPfromQrandom. When the line is put
%in commment this arguement changes at each iterations so that Q is
%randomized. When the line is activated Q is the same for all the
%iterations.
% Lq = [];
% Lq = 1e-1*rand(n);

%% processing
tic

%taking a starting point b0 to initialize the algorithm
bi = PolyFindStartingPoint(A,e_lon,isdiscrete,S);
% bi = bi+i*ones(n,1)*1e-9
Bi = Z - bi*en';
if compute
    %for Nicolas M.
    s(2).X = Bi;
    [stock, Bmat] = stockIt(stock,Bmat,Bi);
    objval = [objval; norm(A-Bi,'fro')];
end

%checking if the problem is real or not, pbisreal is used for the
%computation of the hessian of the barrier function.
if isreal(a)%&& isreal(bi)
    pbisreal = 1;
else
    pbisreal =0;
end

%entering the loop that will iterate to find the best B
h = a-bi;%initial step between A and the starting point
H = h*en';
Dp = speye(n);

try
    while ((norm(H,'fro') > tol)&& (itere < NbItereMax))

        %counting iteration
        itere = itere + 1;

    %---%Find best P and Q
        %same version as the matrix case
        [Dp,Dq,T,invT] = MatFindPQ(Bi,n,e_lon,isdiscrete,verbose);
        %The next line should be actived if Q has to be random
%         [Dp,Dq,T,invT] = MatFindPfromQrandom(Bi,Lq,n,e_lon,1,isdiscrete,verbose);
        P = T*Dp*T';
        iQ = invT'*inv(Dq)*invT;
        s(itere+1).P = P;
        s(itere+1).iQ = iQ;
        P = (P + P')/2;
        iQ = (iQ+iQ')/2;
    
    %---Computing the hessian of the barrier function
        [F] = PolyFindF_old(Z-bi*en',P,iQ,pbisreal,isdiscrete);
        if pbisreal
            am  = a;
            bim = bi;
            sizeF = n;
        else
            am  = [real(a); imag(a)];
            bim = [real(bi);imag(bi)];
            sizeF = 2*n;
        end
        
    %---It is possible to compute the eigenvalues of F so
        %that the problem is easier.
        [U,Df] = eig(F);

        %independent term
        xa = bim-am;
        xa = U'*xa;

    %---%Find the best lambda. The calculus are not depedent on the kind of
    %system (discrete or continuous)
        [la] = PolyFindLambda(xa,Df,eye(sizeF),la,sizeF,tol/1000,verbose);

    %---%Solving the main problem
        h = (eye(sizeF)+la*Df)\xa;
        h = U*h;
%         temp = diag(F);
%         vech = [vech temp];
        if ~pbisreal
            h = complex(h(1:n),h(n+1:2*n));
        end
        H = h*en';

    %---%Updating Bi
        Biprev = Bi;
        bi = bi-h;
        Bi = Z-bi*en';
        

        %         absBi = abs(Bi-Bi')

    %---%saving the evolution of the solution
        if compute
            s(itere+2).X = Bi;
            [stock,Bmat] = stockIt(stock,Bmat,Bi);
            objval = [objval;norm(A-Bi,'fro')];
        end

    %---%print information if necessary
        if verbose
            fprintf('-> Iteration %4i : Distance to the original matrix : %8.6e\n',itere,norm(A-Bi,'fro'))
            fprintf('                            Norm of the direction H : %3.2e\n',norm(H,'fro'))
            fprintf('                            Condition number for P  : %3.2e\n',condest(Dp))
            fprintf('                            Condition number for Q  : %3.2e\n',condest(Dq))
            %fprintf('                    Condition number for T : %3.2e\n',cond(T));
            la %#ok<NOPRT>
            %T %#ok<NOPRT>
            H %#ok<NOPRT>
            Bi %#ok<NOPRT>
        end
    end
 
    toc
    
%% postprocessing
    %returning a polynomial
    bi = poly(Bi);

    %opening figure for the evolution of the solution
    if compute
        close all
        plot_iterates(stock,isdiscrete,e_lon,roc,savf,path);
        save stock.mat stock Bmat objval
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
    mess = err.message;
    id = err.identifier;
    stack = err.stack(1,1);
    err2throw = struct('identifier',id,'stack',stack,'message',mess);
    if compute
        if (strcmp(id,'control:InputArguments')) && (itere>1)
            stock(end,:) = [];
            Bmat(:,end) = [];
            objval(end) = [];
        end
        close all
        plot_iterates(stock,isdiscrete,e_lon,roc,savf,path)
        save stock.mat stock Bmat objval;
        display('Evolution of the Bi saved in stock.mat')
    end
    error(err2throw)
end

