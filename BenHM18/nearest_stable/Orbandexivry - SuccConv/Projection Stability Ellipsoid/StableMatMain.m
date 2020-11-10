function [Bi,e,t]= StableMatMain(A,NbItereMax,timemax,S,isdiscrete,compute,verbose,e_lon,roc,savf,path)
%[Bi]= StableMatMain(A,S,isdiscrete,compute,verbose,e_lon,roc,savf,path)
%StableMain find the closest stable matrix to A in the sense of the
%frobenius norm without explicitly computing the eigenvalues.
%
%   [Bi] = StableMatMain(A,S,isdiscrete,compute,verbose) returns a matrix Bi
%          of the same size than A but with stable eigenvalues. If A has
%          some of its eigenvalues stable, those eigenvalues can be
%          changed. The matrix Bi is computed as the local minimum of an
%          optimization problem.
%
%          The algorithm also tries to keep Bi in the same form than A
%          (diagonal, hermitian, normal,...).
%
%          The algorithm stops when the change between two successive
%          iterations gets smaller than 1e-4(default). If this criterion
%          cannot be satisfied then the algorithm exits with an error.
%
%   Input: A : a square matrix
%          isdiscrete : 1 if A describes a discrete time dynamic system
%                       0 if A describes a continuous time dynamic system
%          compute    : if 1 the algorithm will compute and keep track of
%          the eigenvalues of the system in the file stock.mat. In that
%          file the first line gives the original eigenvalues, the second
%          line gives the eigenvalues of the starting point matrix S ( =
%          -eye(n) by default), the following lines gives the evolution of
%          the eigenvalues. A graph showing the computed values is produced
%          as well. (Blue circles : starting point, red circles : initial
%          eigenvalue, darkness gives the evolution of the iterations).
%          verbose : if 1 print information during the computation
%          S : a starting point matrix whose dimension is equal to A's and
%          whose eigenvalues must be stable.( =-eye(n)+e_lon by default)
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
%  Output: Bi : a stable matrix whose eigenvalues are the closest to A's.
%
%History : 05/01/2011 : Cleaning of the unused code from the file, a backup
%                       (in the bup file dated 31/12/2010) exists which
%                       features a history of the resolution method
%                       'Find_PH'.  +chgt de typo : les fichiers prennent
%                       la particule 'Mat' devant leur nom.
%          18/07/2011 : The method is modified so that solves the problem
%                       in a direct way using Kronecker products. 
%          20/05/2012 : A barrier term is added to enhance (and guarantee)
%                       the convergence.
%          24/08/2012 : Cleaning up of the file after a full backup. In
%                       particular, we don't norm A before using it anymore
%                       since it did not bring anything in precision, it
%                       was difficult to handle (a bug was still
%                       somewhere).
%Last modified : 24/08/2012

%% preprocessing
%adding the Utils directory to matlab path.
% d = pwd;
% smotd = length(d);
% cd ..
% s =pwd;
% smots = length(s);
% d = d(smots+2:smotd);
% cd(d);
% s = strcat(s,'\Utils');
% addpath(s);
% ss = strcat(s,'\plot2svg');
% addpath(ss);
% ss = strcat(s,'\USVD');
% addpath(ss);

cput = cputime; 
%preprocessing
if nargin < 2
    NbItereMax = 100;
end
if nargin < 3
    timemax = +Inf;
end
if nargin < 4
    n = size(A,1); 
    J = (A-A')/2; 
    R = projectPSD(J-A); 
    S = (J-R); 
%S = []; 
end
if nargin < 5
    isdiscrete = 0;
end
if nargin < 6
    compute = 0;
end
if nargin < 7
   verbose = 0;
end
if nargin < 8
    e_lon = 0;
end
if nargin < 9
    roc  = 0;%rate of convergence
    savf = 'none';%type of save
    path = '';
end
if nargin == 0
    MatStableData()
end
if (nargin > 6 && nargin < 9 ) || nargin > 10
    error('-> Wrong number of input arguments')
end


if strcmpi(savf,'none')~= 1
    if ~isdir(path)
        error('The given directory does not exist')
    end
end

if size(A,1)== size(A,2)
    n = size(A,1);
else
    error('-> Input matrix A is not square');
end

%% processing
tic

%parameters
ellips_radius = 0.9999;
beta = 1+1/sqrt(n);%%WARNING : En fonction des problèmes, les valeurs
                          %%de beta 1.01 et 1+1/sqrt(n) apportent des solutions
 %1.01:benchmark          %%différentes. En aucun cas l'une ou l'autre de
 %1+1/sqrt(n) :le reste   %%ces possibilités n'amènent à une solution
                          %%automatiquement optimale. Dans le cas convex,
                          %%1+1/sqrt(n) is optimal pour matrices
tol = n^2*1e-8;%don't put it smaller than 1e-10 otherwise it might be under machine precision in its search of lambda
%tol = 1e-10;
rel_mu_min = n^2*1e-11;%should not be higher than tol
deltatolmu = n^2*1e1;
if tol < 1e-10
    tol = 1e-10;
    warning('Cannot set too small relative accuracy :relative accuracy set to 1e-10'); %#ok<WNTAG>
end

%initialization
itere = 0;
la = 1;

%saving the iterations
if compute
    stock = [];
    Bmat = [];
    [stock, Bmat] = stockIt(stock,Bmat,A);
    VecSqrtObjval = [];
end


%taking a starting point B0 to initialize the algorithm
Bi = MatFindStartingPoint(A,e_lon,isdiscrete,S);


%Storage of the values
if compute
    [stock, Bmat] = stockIt(stock,Bmat,Bi);
    VecSqrtObjval = [VecSqrtObjval; norm(A-Bi,'fro')];
end

%Central path parameter
lmin = max(real(eig(Bi)));
mu  = min(1/sqrt(n*lmin^2),sqrt(n*lmin^2));

% [Dp,Dq,T,invT] = MatFindPQ(Bi,n,e_lon,isdiscrete,verbose);
% iQ = invT'*inv(Dq)*invT;
% P = T*Dp*T';
% iQP = 2*invT'*inv(Dq)*Dp*T';
% bsx2 = iQ*(iQP*P+P*iQP')*iQP;
% mu = -real(trace(bsx2'*(Bi-A)))/(real(trace(bsx2'*iQP)))

%entering the loop that will iterate to find the best B
Y = Bi;
logdetY = -log(real(det(Y)));
sqrtObjval = norm(A-Bi,'fro');
H = sqrtObjval;
DiffNormObjval = 0.5*sqrtObjval^2;%initial step between A and the starting point
DiffObjFunction = DiffNormObjval; %negating the logdet contribution ensures that it goes at least once in the while loop)
mumin = min(abs(rel_mu_min*DiffNormObjval),abs(rel_mu_min*DiffNormObjval/logdetY));

try
    %while (((DiffObjFunction > tol)|| (mu~=mumin) )&& (itere < NbItereMax))
    while (itere < NbItereMax) && cputime-cput <= timemax 
        %         if norm(H,'fro') < 1e-1%for debug use
        %             display('under first tol')
        %         end

        %counting iteration
        itere = itere + 1;

    %---%Find best P and Q
        try
            [Dp,Dq,T,invT] = MatFindPQ(Bi,n,e_lon,isdiscrete,verbose);
        catch E 
            disp('The algorithm has run into numerical problems:');
            disp(E); 
            return; 
        end
        
    %---%Construct the system in kronecker form for the direct method
        [F,w] = MatDirectFindF(A,Bi,Dp,Dq,T,invT,n,mu,isdiscrete,verbose);
        sizeF = 2*n^2;
        
        
    %---It is possible to compute the eigenvalues of F so
        %that the problem is easier.
        [U,Df] = eig(F);
        Dg = eye(sizeF)+mu*Df;

        %independent term
        ax = U'*w;
        
        
    %---%Computing Lambda
        [la] = PolyFindLambda(ax,Df,Dg,la,sizeF,tol*1e-3,ellips_radius,verbose);

    %---%FindPH_approx : Does not try to solve directly the equation of the
        %projection on the ellipsoid but rather tries to find an
        %approximation to it. It norms 'G', the derivative of the distance
        %to A. Since the two derivatives need only to be equal in
        %direction, not in norm. This correspond at looking for a point on
        %an ellipsoid of unit radius. It also features a stopping criterion
        %based on an inner ellipsoid of smaller radius.
        %The stopping criterion based on the inner ellipsoid was achieved
        %only for the discrete case. The file has also been cleaned of
        %debug operations. 
        %Creation :23/09/2010 
        %Modifications : 18/07/2011 : Use of a direct method for creating
        %examples for the article
        %Last update:18/07/2011
%         [H11] = MatFindPH_approx_updated(Dp,Dq,Bi,A,T,invT,n,la,isdiscrete,e_lon,gamma,verbose);
%         [H12] = MatFindPH_approx_conditioned(Dp,Dq,Bi,A,T,invT,n,la,isdiscrete,e_lon,gamma,verbose);
        uh = ax./(diag(Dg)+la*diag(Df));%(eye(sizeF)+la*Df)\ax; same result as the other
%         uhs = axs./(ones(sizeF,1)+las*diag(Dfs));
        h = U*uh;
%         hs = Us*uhs;
        h = complex(h(1:n^2),h(n^2+1:2*n^2));
%         hs = complex(hs(1:n^2),hs(n^2+1:2*n^2));
        H = reshape(h,n,n);
%         Hs = reshape(hs,n,n);
        
%                 diffH = H10-H9
%                 H10/H9
%                 norm(diffH,'fro')/norm(H,'fro')
        %         DiffH2H1 = H5-H4 %#ok<NASGU>
        %         h_1 = (H3+H3')/2;
        %         h_2 = (H4+H4')/2;
        %         diffh2h1 = h_2-h_1
        %         a_1 = (H3-H3')/2;
        %         a_2 = (H4-H4')/2;
        %         diffa2a1 = a_2-a_1

    %---%Updating Bi
        Bi = Bi+H;        
        
    %---%Computing the update of the stopping criterion
        sqrtObjvalprev = sqrtObjval;
        sqrtObjval = norm(A-Bi,'fro');
        DiffNormObjval = 0.5*(sqrtObjvalprev^2-sqrtObjval^2);
        logdetYprev = logdetY;
        P = T*Dp*T';
        if ~isdiscrete
            Y = Bi*P;
            Y = -(Y+Y');
            %Q = Y;
        else
            Y = [P 2*P*Bi; zeros(size(Bi)) P];
            Y = (Y+Y')/2;
            %Q = P-Bi*P*Bi';
        end
        logdetY = -log(real(det(Y)));    
        DiffObjFunction = abs(DiffNormObjval + mu*(logdetYprev-logdetY));
        
    %--Update of mu
        if DiffObjFunction < tol*deltatolmu
            mu = max(mu/beta,mumin);
        end

    %---%saving the evolution of the solution
        if compute
            [stock,Bmat] = stockIt(stock,Bmat,Bi);
            VecSqrtObjval = [VecSqrtObjval;sqrtObjval]; %#ok<AGROW>
        end

    %---%print information if necessary
        if verbose
            fprintf('\n-> Iteration %4i : Distance to the original matrix : %8.6e\n',itere,sqrtObjval)
            fprintf('                     Improvement of the obj. val.   : %3.2e\n',DiffObjFunction)
            fprintf('                     Norm of the direction H        : %3.2e\n',norm(H,'fro'))
            fprintf('                     Value of the barrier param.    : %3.2e\n',mu)
%             fprintf('cond P : %d\n',cond(P)) 
%             fprintf('cond Q : %d\n',cond(Q)) 
%             H %#ok<NOPRT>
%             Bi %#ok<NOPRT>
        end
        e(itere) = norm(A-Bi,'fro'); 
        t(itere) = cputime-cput; 
        if mod(itere,10) == 0
            fprintf('%2.0f:%2.3f - ',itere,e(itere)); 
        end
        if mod(itere,100) == 0
            fprintf('\n');       
        end
    end
    fprintf('\n');    
    
    %opening figure for the evolution of the solution
    if compute
        close all
        plot_iterates(stock,isdiscrete,e_lon,roc,savf,path);
        VecSqrtObjval = [min(VecSqrtObjval);VecSqrtObjval];
        save stock.mat stock Bmat VecSqrtObjval
    end
    
    %convergence plot
    if compute && verbose
        %capting the Obective function's values
        noAX = VecSqrtObjval(2:end).^2;%VecObjVal(1) contains the distance to the starting point

        %first plot : Evolution of the Objective function
        figure;
        hold on;
        plot(noAX,'g')
        plot(noAX,'go')
        title('Norm A-X');
        hold off;

        %second plot : plot of the convergence
        %second plot : plot of the convergence
        figure;
        nk = length(noAX);
        k = nk;
        %     loglog(1:k,noAX(1:k)-noAX(nk),'r')
        loglog(noAX(1:k-100)-noAX(k),'r')
        hold on;
        t = logspace(0,log10(itere));
        lt = length(t);
        loglog(t,t.^(-1))
        th = t(1:round(lt/2));
        loglog(th,th.^(-2))
        hold off
    end

    %displaying time
    time = toc;
%     if verbose
%        time
%     end

catch
    toc
    %opening figure for the evolution of the solution
%     err = lasterror;
%     msgid = err.identifier;
    if compute
%         if (strcmp(msgid,'control:InputArguments')) && (itere>1)
%             stock(end,:) = [];
%             Bmat(:,end) = [];
%             objval(end) = [];
%         end

        %close all
        plot_iterates(stock,isdiscrete,e_lon,roc,savf,path);
        VecSqrtObjval = [min(VecSqrtObjval);VecSqrtObjval];
        save stock.mat stock Bmat VecSqrtObjval
        display('Evolution of the Bi saved in stock.mat')
    end
    if compute && verbose
        %capting the Obective function's values
        noAX = VecSqrtObjval(2:end).^2;%VecObjVal(1) contains the distance to the starting point

        %first plot : Evolution of the Objective function
        figure;
        hold on;
        plot(noAX,'g')
        plot(noAX,'go')
        title('Norm A-X');
        hold off;

        %second plot : plot of the convergence
        figure;
        nk = length(noAX);
        k = nk;
        %     loglog(1:k,noAX(1:k)-noAX(nk),'r')
        loglog(noAX(1:k-100)-noAX(k),'r')
        hold on;
        t = logspace(0,log10(itere));
        lt = length(t);
        loglog(t,t.^(-1))
        th = t(1:round(lt/2));
        loglog(th,th.^(-2))        
        hold off
        %second plot : plot of the convergence

    end
    rethrow(lasterror)
end

