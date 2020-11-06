function RUNME()
%%% summary: GAMM / PAMM 2017 companion code
%%% project: Fast Low-Rank Empirical Cross Gramians
%%% authors: Christian Himpe ( 0000-0003-2194-6754 )
%%% license: 2-Clause BSD (2017)
%$
    if(exist('emgr')~=2)
        error('emgr not found! Get emgr at: http://gramian.de');
    else
        global ODE;
        ODE = [];
        fprintf('emgr (version: %1.1f)\n',emgr('version'));
    end

%% SETUP

    rand('seed',1009);

    K = 1024;
    del = rand(K,1)*0.001;
    omg = rand(K,1)*100.0;
    a = cell(K);
    for k=1:K
        a{k} = sparse([-2.0*del(k)*omg(k),-omg(k);omg(k),0]);
    end
    A = blkdiag(a{:});

    N = 2*K;
    M = 1;
    Q = M;

    B = ones(N,1);
    C = ones(1,N).*( kron(1.0./omg',[0,1.0]) + kron(ones(1,K),[1.0,0]) );

    T = [0.001,1.0];
    L = floor(T(2)/T(1)) + 1;
    U = @(t) ones(M,1)*(t<=T(1))/T(1);
    X = zeros(N,1);

    LIN = @(x,u,p,t) A*x + B*u;
    OUT = @(x,u,p,t) C*x;

%% FULL ORDER
    Y = ODE(LIN,OUT,T,X,U,0);
   %figure; plot(0:T(1):T(2),Y); return;
    n1 = norm(Y(:),1);
    n2 = norm(Y(:),2);
    n8 = norm(Y(:),Inf);

%% OFFLINE
    global DWX;
    DWX = [];
    pre = 1e-13;
%
    o1 = tic;
    WX = emgr(LIN,OUT,[M,N,Q],T,'x');
    OFFLINE_FULL_WX = toc(o1)
    o1 = tic;
    [U1,D1] = hapod({WX},pre,'none');
    OFFLINE_FULL_SV = toc(o1)
    FULL_BASE_SIZE= size(U1,2)
%
    o2 = tic;
    w = 64;
    n = ceil(N/w);
    wx = cell(1,n);
    for k=1:n
        wx{k} = emgr(LIN,OUT,[M,N,Q],T,'x',0,[0,0,0,0,0,0,0,0,0,0,w,k]);
    end;
    OFFLINE_DIST_WX = toc(o2)
    o2 = tic;
    [U2,D2] = hapod(wx,pre,'dist',0.75);
    OFFLINE_DIST_SV = toc(o2)
    DIST_BASE_SIZE = size(U2,2)
    RESIDUAL = norm(WX-cell2mat(wx),2)
%

%% EVALUATION
    m = min(size(U1,2),N-1);
    for n=1:m
        uu = U1(:,1:n);
        a = uu'*A*uu;
        b = uu'*B;
        c = C*uu;
        x = uu'*X;
        lin = @(x,u,p,t) a*x + b*u;
        out = @(x,u,p,t) c*x;
        y = ODE(lin,out,T,x,U,0);
        l1(n) = norm(Y(:)-y(:),1)/n1;
        l2(n) = norm(Y(:)-y(:),2)/n2;
        l8(n) = norm(Y(:)-y(:),Inf)/n8;
    end;

    twx = [OFFLINE_FULL_WX,OFFLINE_DIST_WX];
    tsv = [OFFLINE_FULL_SV,OFFLINE_DIST_SV];

    save('runme.mat','l1','l2','l8','m','twx','tsv');
end

