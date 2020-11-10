function [S,U,V,Q] = USVD(A,B,varargin)
% [S,U,V,Q] = USVD(A,B,'option',value)computes the Unitary Singular Value
% Decomposition of C = A'*B by applying transformation to the elements of
% the product separatly.
% Input  : A and B two n by n and non-singular matrices. The matrix A is
%          lower triangular and B is upper triangular.
%          optional input : 'tol' to set the accuracy of the method
%                           'verbose' set to 1 or 0
% Output : S and unitary U,V and Q such that U'*A'*B*V = S and U'*A'*Q' and
%          Q*B*V are upper triangular.

[m p] = size(B);
[q n] = size(A);

%checking square and same size
if (m ~= p)||(q ~= n)||(m ~= n)
    error('Input matrices A and B are not square and of the same size')
end

%parameters
tol = 1e-16;
verbose = 0;

m = length(varargin);
if m > 4
    error('Number of input is to big')
elseif m == 2
    if strcmpi(varargin{1},'tol')
        tol = varargin{2};
    end
    if strcmpi(varargin{1},'verbose')
        verbose = varargin{2};
    end
elseif m == 4
    if strcmpi(varargin{1},'tol')
        tol = varargin{2};
        verbose = varargin{4};
    elseif strcmpi(varargin{1},'verbose')
        verbose = varargin{2};
        tol = varargin{4};
    else
        error('Unknown input name')
    end
end

%initialization
k = 1;
Uk = eye(n);
Vk = eye(n);
Qk = eye(n);
saveTrace = 0;
achieved_precision = tol +1;
wquit = 0;


while (achieved_precision > tol) && ~wquit

    %sweep when upper triangular
    [A,B,Up,Vp,Qp] = sweep(A,B,n);

    %sweep when lower triangular (this is equivalent to make the sweep to
    %(B,A) = B'*A
    [B,A,Vl,Ul,Ql] = sweep(B,A,n);

    %computing transformations
    Uk = Uk*Up*Ul;
    Vk = Vk*Vp*Vl;
    Qk = Ql*Qp*Qk;

    Sk = A'*B;

    %records
    saveTracePrev = saveTrace;
    saveTrace = trace(Sk);
    achieved_precision = abs(saveTrace-saveTracePrev);

    k = k + 1;
    if (k >= n)&&(k>6)
        if verbose
            warning('Process seems to not converge in USVD.m')
        end
        wquit = 1;
    end
end

%termination
S = Sk;
U = Uk;
V = Vk;
Q = Qk;

