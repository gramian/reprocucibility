function RUNME(s)
%%% summary: Companion Code to "On Reduced Input-Output Dynamic Mode Decomposition"
%%% authors: Christian Himpe ( 0000-0003-2194-6754 ),
%%%          Tim Mitchell ( 0000-0002-8426-0242 )
%%% license: 2-Clause BSD (2017)
%
% Usage:
%
%  RUNME(0) : run reduced ioDMD examples without stabilization
%
%  RUNME(1) : run reduced ioDMD examples with stabilization
%$
    addpath('GRANSO');

    if(nargin==0 || isempty(s)), s = 1; end;

    % Setup System
    randn('seed',1009);                         % seed random numbers
    nx = 1000;                                  % state space dimension
    nt = nx;                                    % number of time steps
    T = 1.0;                                    % time horizon
    dt = T./nt;                                 % time step width
    [A,B,C] = make_1D_linear_transport(nx);     % assemble semi discrete model

    % Generate target data
    X0 = zeros(nx,1);                           % initial state
    U0 = @(t) exp(((t-0.1).^2)./(-0.001));      % input function
    [Ua,Xa,Ya] = rk1i(A,B,C,0,[dt,T],X0,U0);    % compute input. state, output functions
    n0 = norm(Ya(:),2);                         % compute discrete L2 norm

    % Generate sampling data (Random input)
    [Ub,Xb,Yb] = rk1i(A,B,C,0,[dt,T],zeros(nx,1),@(t) randn(1));

    % Generate sampling data (Random initial, intermediary function for cross excitation)
    [Uc,Xc,Yc] = rk1i(A,B,C,0,[dt,T],randn(nx,1),@(t) 0);

    % Generate sampling data (Cross excitation via random initial state)
    [Ud,Xd,Yd] = rk1i(A,B,C,0,[dt,T],zeros(nx,1),@(t) Yc(:,1+floor(t/dt)));

    % Generate sampling data (Step input)
    [Ue,Xe,Ye] = rk1i(A,B,C,0,[dt,T],zeros(nx,1),@(t) 1);

    % Generate sampling data (Shifted initial, intermediary function for cross excitation)
    [Uf,Xf,Yf] = rk1i(A,B,C,0,[dt,T],ones(nx,1),@(t) 0);

    % Generate sampling data (Cross excitation, via shifted initial state)
    [Ug,Xg,Yg] = rk1i(A,B,C,0,[dt,T],zeros(nx,1),@(t) Yf(:,1+floor(t/dt)));

    l1 = 0;    % No regularization
    l2 = 1e-5; % Regularization coefficient

    K = 8;
    oe = zeros(K,10); % relative output error
    ot = zeros(K,10); % stabilization runtime
    os = zeros(K,5);  % order of reduced models
    st = zeros(K,10); % stabilzation test

    data = {};

    for k=1:K

        fprintf(2,'\nRunning %d of %d ...\n',k,K);

        r = 10^-k; % prescribed projection error

        xa = pod(Xa,r)'*Xa; os(k,1) = size(xa,1); % compress via POD

        % unregularized ioDMD using target data
        [A0,B0,C0,D0] = iodmd(Ua,xa,Ya,l1);
        if(s), z = tic; [A0,B0,C0,D0,data{end+1}] = stabilize(A0,B0,C0,D0,Ua,xa,Ya); ot(k,1) = toc(z); end;
        Y = dint(A0,B0,C0,D0,[dt,T],zeros(size(A0,2),1),U0);
        oe(k,1) = norm(Ya(:)-Y(:),2)/n0;
        st(k,1) = all(abs(eig(A0))<1);

        % regularized ioDMD using target data
        [A0,B0,C0,D0] = iodmd(Ua,xa,Ya,l2);
        if(s), z = tic; [A0,B0,C0,D0,data{end+1}] = stabilize(A0,B0,C0,D0,Ua,xa,Ya); ot(k,2) = toc(z); end;
        Y = dint(A0,B0,C0,D0,[dt,T],zeros(size(A0,2),1),U0);
        oe(k,2) = norm(Ya(:)-Y(:),2)/n0;
        st(k,2) = all(abs(eig(A0))<1);

        xb = pod(Xb,r)'*Xb; os(k,2) = size(xb,1); % compress via POD

        % unregularized ioDMD using random noise input
        [A0,B0,C0,D0] = iodmd(Ub,xb,Yb,l1);
        if(s), z = tic; [A0,B0,C0,D0,data{end+1}] = stabilize(A0,B0,C0,D0,Ub,xb,Yb); ot(k,3) = toc(z); end;
        Y = dint(A0,B0,C0,D0,[dt,T],zeros(size(A0,2),1),U0);
        oe(k,3) = norm(Ya(:)-Y(:),2)/n0;
        st(k,3) = all(abs(eig(A0))<1);

        % regularized ioDMD using random noise input
        [A0,B0,C0,D0] = iodmd(Ub,xb,Yb,l2);
        if(s), z = tic; [A0,B0,C0,D0,data{end+1}] = stabilize(A0,B0,C0,D0,Ub,xb,Yb); ot(k,4) = toc(z); end;
        Y = dint(A0,B0,C0,D0,[dt,T],zeros(size(A0,2),1),U0);
        oe(k,4) = norm(Ya(:)-Y(:),2)/n0;
        st(k,4) = all(abs(eig(A0))<1);

        xd = pod(Xd,r)'*Xd; os(k,3) = size(xd,1); % compress via POD

        % unregularized ioDMD using cross excitation w. random initial state
        [A0,B0,C0,D0] = iodmd(Ud,xd,Yd,l1);
        if(s), z = tic; [A0,B0,C0,D0,data{end+1}] = stabilize(A0,B0,C0,D0,Ud,xd,Yd); ot(k,5) = toc(z); end;
        Y = dint(A0,B0,C0,D0,[dt,T],zeros(size(A0,2),1),U0);
        oe(k,5) = norm(Ya(:)-Y(:),2)/n0;
        st(k,5) = all(abs(eig(A0))<1);

        % regularized ioDMD using cross excitation w. random initial state
        [A0,B0,C0,D0] = iodmd(Ud,xd,Yd,l2);
        if(s), z = tic; [A0,B0,C0,D0,data{end+1}] = stabilize(A0,B0,C0,D0,Ud,xd,Yd); ot(k,6) = toc(z); end;
        Y = dint(A0,B0,C0,D0,[dt,T],zeros(size(A0,2),1),U0);
        oe(k,6) = norm(Ya(:)-Y(:),2)/n0;
        st(k,6) = all(abs(eig(A0))<1);

        xe = pod(Xe,r)'*Xe; os(k,4) = size(xe,1); % compress via POD

        % unregularized ioDMD using constant excitation
        [A0,B0,C0,D0] = iodmd(Ue,xe,Ye,l1);
        if(s), z = tic; [A0,B0,C0,D0,data{end+1}] = stabilize(A0,B0,C0,D0,Ue,xe,Ye); ot(k,7) = toc(z); end;
        Y = dint(A0,B0,C0,D0,[dt,T],zeros(size(A0,2),1),U0);
        oe(k,7) = norm(Ya(:)-Y(:),2)/n0;
        st(k,7) = all(abs(eig(A0))<1);

        % regularized ioDMD using constant excitation
        [A0,B0,C0,D0] = iodmd(Ue,xe,Ye,l2);
        if(s), z = tic; [A0,B0,C0,D0,data{end+1}] = stabilize(A0,B0,C0,D0,Ue,xe,Ye); ot(k,8) = toc(z); end;
        Y = dint(A0,B0,C0,D0,[dt,T],zeros(size(A0,2),1),U0);
        oe(k,8) = norm(Ya(:)-Y(:),2)/n0;
        st(k,8) = all(abs(eig(A0))<1);

        xg = pod(Xg,r)'*Xg; os(k,5) = size(xg,1); % compress via POD

        % unregularized ioDMD using cross excitation w. shifted initial state
        [A0,B0,C0,D0] = iodmd(Ug,xg,Yg,l1);
        if(s), z = tic; [A0,B0,C0,D0,data{end+1}] = stabilize(A0,B0,C0,D0,Ug,xg,Yg); ot(k,9) = toc(z); end;
        Y = dint(A0,B0,C0,D0,[dt,T],zeros(size(A0,2),1),U0);
        oe(k,9) = norm(Ya(:)-Y(:),2)/n0;
        st(k,9) = all(abs(eig(A0))<1);

        % regularized ioDMD using cross excitation w. shifted initial state
        [A0,B0,C0,D0] = iodmd(Ug,xg,Yg,l2);
        if(s), z = tic; [A0,B0,C0,D0,data{end+1}] = stabilize(A0,B0,C0,D0,Ug,xg,Yg); ot(k,10) = toc(z); end;
        Y = dint(A0,B0,C0,D0,[dt,T],zeros(size(A0,2),1),U0);
        oe(k,10) = norm(Ya(:)-Y(:),2)/n0;
        st(k,10) = all(abs(eig(A0))<1);
    end

    save('data.mat','data');

    fprintf('\n');

    % Plot runtime of stabilization for different state compression levels
    if(s)
        figure;
        semilogy(1:K,ot(:,1),'r-*','linewidth',3,'markersize',5);
        hold on;
        semilogy(1:K,ot(:,3),'g-*','linewidth',3,'markersize',5);
        semilogy(1:K,ot(:,5),'b-*','linewidth',3,'markersize',5);
        semilogy(1:K,ot(:,7),'k-.*','linewidth',3,'markersize',5);
        semilogy(1:K,ot(:,9),'m--*','linewidth',3,'markersize',5);
        hold off;
        legend('Target','PE (Gauss)','CE (Gauss)','PE (Step)',' CE (Shift)','location','northwest');
        xlim([1,K]);
        ylim([0.05,100]);
        set(gca,'XTickLabel',10.^-(2:2:K));
        xlabel('Projection Error');
        ylabel('Stabilization Runtime');
        pbaspect([1,1,1]);
        set([gca; findall(gca,'Type','text')],'FontSize',16);
        print('-depsc',['prop',int2str(1+(s~=0)),'b.eps']);
    else
        st
    end

    % Plot reduced order for different state compression levels
    figure;
    plot(1:K,os(:,1),'r-*','linewidth',3,'markersize',5);
    hold on;
    plot(1:K,os(:,2),'g-*','linewidth',3,'markersize',5);
    plot(1:K,os(:,3),'b-*','linewidth',3,'markersize',5);
    plot(1:K,os(:,4),'k-*','linewidth',3,'markersize',5);
    plot(1:K,os(:,5),'m--*','linewidth',3,'markersize',5);
    hold off;
    legend('Target','PE (Gauss)','CE (Gauss)','PE (Step)',' CE (Shift)','location','northwest');
    xlim([1,K]);
    set(gca,'XTickLabel',10.^-(2:2:K));
    xlabel('Projection Error');
    ylabel('System Dimension');
    pbaspect([1,1,1]);
    set([gca; findall(gca,'Type','text')],'FontSize',16);
    print('-depsc',['prop',int2str(1+(s~=0)),'a.eps']);

    % Plot relative discrete L2 output error for different state compression levels (unregularized)
    figure;
    semilogy(1:K,oe(:,1),'-r*','linewidth',3,'markersize',5);
    hold on;
    semilogy(1:K,oe(:,3),'-g*','linewidth',3,'markersize',5);
    semilogy(1:K,oe(:,5),'-b*','linewidth',3,'markersize',5);
    hold off;
    set(gca,'XTickLabel',10.^-(2:2:K));
    xlim([1,K]);
    ylim([1e-3,2]);
    xlabel('Projection Error');
    ylabel('Output Error');
    pbaspect([1,1,1]);
    set([gca; findall(gca,'Type','text')],'FontSize',16);
    legend('Target','PE','CE','location','southeast');
    print('-depsc',['numex',int2str(1+(s~=0)),'a.eps']);

    % Plot relative discrete L2 output error for different state compression levels (reguarized)
    figure;
    semilogy(1:K,oe(:,2),'-r*','linewidth',3,'markersize',5);
    hold on;
    semilogy(1:K,oe(:,4),'-g*','linewidth',3,'markersize',5);
    semilogy(1:K,oe(:,6),'-b*','linewidth',3,'markersize',5);
    hold off;
    set(gca,'XTickLabel',10.^-(2:2:K));
    xlim([1,K]);
    ylim([1e-3,2]);
    xlabel('Projection Error');
    ylabel('Output Error');
    pbaspect([1,1,1]);
    set([gca; findall(gca,'Type','text')],'FontSize',16);
    legend('Target','PE','CE','location','southeast');
    print('-depsc',['numex',int2str(1+(s~=0)),'c.eps']);

    % Plot relative discrete L2 output error for different state compression levels (unregularized)
    figure;
    semilogy(1:K,oe(:,1),'-r*','linewidth',3,'markersize',5);
    hold on;
    semilogy(1:K,oe(:,7),'-g*','linewidth',3,'markersize',5);
    semilogy(1:K,oe(:,9),'--b*','linewidth',3,'markersize',5);
    hold off;
    set(gca,'XTickLabel',10.^-(2:2:K));
    xlim([1,K]);
    ylim([1e-12,2]);
    xlabel('Projection Error');
    ylabel('Output Error');
    pbaspect([1,1,1]);
    set([gca; findall(gca,'Type','text')],'FontSize',16);
    legend('Target','PE','CE','location','southwest');
    print('-depsc',['numex',int2str(1+(s~=0)),'b.eps']);

    % Plot relative discrete L2 output error for different state compression levels (regularized)
    figure;
    semilogy(1:K,oe(:,2),'-r*','linewidth',3,'markersize',5);
    hold on;
    semilogy(1:K,oe(:,8),'-g*','linewidth',3,'markersize',5);
    semilogy(1:K,oe(:,10),'--b*','linewidth',3,'markersize',5);
    hold off;
    set(gca,'XTickLabel',10.^-(2:2:K));
    xlim([1,K]);
    ylim([1e-12,2]);
    xlabel('Projection Error');
    ylabel('Output Error');
    pbaspect([1,1,1]);
    set([gca; findall(gca,'Type','text')],'FontSize',16);
    legend('Target','PE','CE','location','southwest');
    print('-depsc',['numex',int2str(1+(s~=0)),'d.eps']);
end

% Linear Transport Assembly
function [A,B,C] = make_1D_linear_transport(nx)

    a = 1.3;
    A = (a*nx)*spdiags([ones(nx,1),-ones(nx,1)],[-1,0],nx,nx);
    A(1,1) =1e1 -nx;
    B = [nx*1.0;sparse(nx-1,1)];
    C = [sparse(1,nx-1),1.0];
end

% Proper Orthogonal Decomposition
function M = pod(data,bound)

    [U,D,V] = svd(data,'econ');

    D = diag(D);
    d = flipud(cumsum(flipud(D.^2)));

    K = find(d <= bound^2,1);
    K = K + (K==1);
    if(isempty(K)), K = min(size(data)) + 1; end;
    M = U(:,1:K-1);
end

% Input-Output Dynamic Mode Decomposition
function [A,B,C,D] = iodmd(U,X,Y,l)

    m = size(U,1);
    n = size(X,1);
    q = size(Y,1);

    [u,d,v] = svd([X(:,1:end-1);U(:,1:end-1)],'econ');

    s = length(find(diag(d)>l));

    u = u(:,1:s);
    v = v(:,1:s);
    di = diag(1.0./diag(d(1:s,1:s)));

    A = X(:,2:end) * v * di * u(1:n,:)';
    B = X(:,2:end) * v * di * u(n+1:n+m,:)';
    C = Y(:,1:end-1) * v * di * u(1:n,:)';
    D = Y(:,1:end-1) * v * di * u(n+1:n+q,:)';
end

% Optimization-Based Stabilization
function [A,B,C,D,data] = stabilize(A,B,C,D,U,X,Y)

    maxit = 1000;
    H = [A,B;C,D];
    prob = makeProblem(U(:,1:end-1),X(:,1:end-1),X(:,2:end),Y(:,1:end-1));
    [halt_log_fn, get_log_fn] = makeHaltLogFunctions(maxit);
    opts = gransoOptions(prob.nvar);
    opts.x0 = H(:);
    opts.maxit = maxit;
    opts.print_level = 0;
    opts.halt_log_fn = halt_log_fn;
    opts.opt_tol = 1e-8;
    %opts.limited_mem_size = 10; % Uncomment to use Limited-Memory-BFGS (not recommended)
    [soln,t] = runGranso(prob,opts);
    soln = rmfield(soln,'H_final'); % This is critical
    data.soln = soln;
    data.log = get_log_fn();
    data.log.t = t;
    [A,B,C,D] = prob.xToSubMatrices(soln.most_feasible.x);
end

% Discrete Integrator
function y = dint(A,B,C,D,t,x0,u)

    h = t(1);
    K = floor(t(2)/h) + 1;

    y = zeros(size(C,1),K);
    y(:,1) = C*x0 + D*u(0);

    xk = x0;

    for k = 2:K
        tk = (k - 1.5) * h;
        uk = u(tk);
        xk = A*xk + B*uk;
        y(:,k) = C*xk + D*uk;
    end;
end

% Implicit First-Order Runge-Kutta
function [U,X,Y] = rk1i(A,B,C,D,t,x,u)

    h = t(1);
    K = floor(t(2)/h) + 1;

    U = u(0);
    U(end,K) = 0;
    X = x;
    X(end,K) = 0;
    Y = C*x;
    Y(end,K) = 0;

    [la,ua,pa] = lu(speye(size(A,1)) - h*A,'vector');

    xk = x;

    for k=2:K
        tk = (k - 1.5) * h;
        uk = u(tk);
        xk = ua\(la\(xk(pa) + h*B*u(tk) ) );
        U(:,k) = uk;
        X(:,k) = xk;
        Y(:,k) = C*xk;
    end
end
