function RUNME()
%%% summary: runme (emgr - the empirical gramian framework companion code)
%%% authors: Christian Himpe ( 0000-0003-2194-6754 )
%%% license: 2-Clause BSD (2016--2018)
%$

    rand('seed',1009);
    randn('seed',1009);

%% NumEx 1

    sys.M = 4;									% Number of inputs
    sys.N = 256;								% Number of states
    sys.Q = 4;									% Number of outputs
    sys.h = 0.01;								% Time step
    sys.T = 1.0;								% Time horizon

    a = 1e-1;
    b = 1e+1;
    WX = -diag( a*((b/a).^rand(sys.N,1)) );					% Balanced cross gramian
    B = randn(sys.N,sys.M);							% Balanced input and output matrix
    A = sylvester(WX,WX,B*B') - sqrt(eps)*speye(sys.N);				% Balanced system matrix
    Q = orth(randn(sys.N,sys.N));						% Unbalancing transformation
    A = Q'*A*Q;									% Unbalanced system matrix
    B = Q'*B;									% Unbalanced input matrix
    C = B';									% Unbalanced output matrix

    sys.f = @(x,u,p,t) A*x + B*u;						% Vector field
    sys.F = @(x,u,p,t) A'*x + C'*u;						% Adjoint vector field
    sys.g = @(x,u,p,t) C*x;							% Output functional

    V = randn(sys.M,floor(sys.T/sys.h) + 1);
    sys.vt = @(t) V(:,floor(t/sys.h)+1);

    curios(sys,'state-reduction','linear-direct-truncation');
    curios(sys,'system-index','hinf-bound',{'hold','cached'});

%% NumEx 2

    sys.M = 1;									% Number of inputs
    sys.N = 256;								% Number of states
    sys.Q = 1;									% Number of outputs
    sys.h = 1.0./sys.N;								% Time step
    sys.T = 1.0;								% Time horizon

    A = sys.N * spdiags(ones(sys.N,1)*[1,-1],[-1,0],sys.N,sys.N);		% System matrix
    A(1,1) = 0;
    B = sparse(1,1,sys.N,sys.N,1);						% Input matrix
    C = sparse(1,sys.N,1.0,1,sys.N);						% Output matrix

    sys.f = @(x,u,p,t) p*A*x + B*u;						% Vector field
    sys.F = @(x,u,p,t) p*A'*x + C'*u;						% Adjoint vector field
    sys.g = @(x,u,p,t) C*x;							% Output functional

    sys.ut = Inf;
    sys.vt = @(t) exp(((t-0.1).^2)./(-0.001));					% Input function
    sys.p = [1.0,1.5];								% Training transport velocities
    sys.q = 0.5*rand(1,10) + 1.0;						% Test transport velocities

    curios(sys,'state-reduction','linear-direct-truncation');

%% NumEx 3

    sys.M = 1;									% Number of inputs
    sys.N = 256;								% Number of states
    sys.Q = 4;									% Number of outputs
    sys.h = 0.01;								% Time step
    sys.T = 1.0;								% Time horizon

    A = -gallery('lehmer',sys.N);						% System matrix
    B = 0.5 + 0.5*cos(2.0*pi/sys.N*[1:sys.N]');					% Input matrix
    C = kron(eye(sys.Q),ones(1,sys.N/sys.Q)); 					% Output matrix

    sys.f = @(x,u,p,t) A*tanh(p.*x) + B*u;					% Vector field
    sys.g = @(x,u,p,t) C*x;							% Output functional
    sys.p = ones(sys.N,1) * [0.5,1.0];						% Training parameter range
    sys.q = 0.5 + 0.5 * rand(sys.N,10);						% Test parameter

    curios(sys,'combined-reduction','minimality-based',{'active','linpar','nonsym'});
end
