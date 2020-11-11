function numex(experiment)
%%% project: Cross-Gramian-Based Dominant Subspaces Companion Code
%%% version: 1.3 ( 2019-08-15 )
%%% authors: C. Himpe ( 0000-0003-2194-6754 )
%%% license: BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%%% summary: Dominant subspaces and balanced truncation for two benchmarks

    global ODE;

    switch(experiment)

        case 2.1 % Convection Benchmark 0

            load('convection0.mat');    % Loads A, B, C, E
            dt = 0.04; 			% time-step width
            Tf = 4.0;  			% time horizon
            epsilon = 10.0.^(-[2:8]);	% Projection errors
            ut = @(t) (t<=dt);		% test input

        case 2.2 % Convection Benchmark 1

            load('convection1.mat');    % Loads A, B, C, E
            dt = 0.04; 			% time-step width
            Tf = 4.0;  			% time horizon
            epsilon = 10.0.^(-[2:8]);	% Projection errors
            ut = @(t) (t<=dt);		% test input

        otherwise % FOM Benchmark

            load('fom.mat');		% Loads A, B, C, E
            E = speye(1006);		% Unit mass matrix
            dt = 0.001; 		% time-step width
            Tf = 1.0;  			% time horizon
            epsilon = 10.0.^(-[3:12]);	% Projection errors
            ut = @(t) (t<=dt);		% test input
    end

    M = size(B,2);								% Number of inputs
    N = size(A,2);								% Number of states
    Q = size(C,1);								% Number of outputs

    F = @(x,u,p,t) 0;								% Vector field
    G = @(x,u,p,t) 1;								% Adjoint vector field

    g = @(x,u,p,t) C*x;								% Output function

    X0 = zeros(N,1);								% Initial condition

    Y = rk1i(A,B,C,E,[dt,Tf],X0,ut,0,0);					% Reference solution
    n2 = norm(Y(:),2);								% L2 norm of reference solution

    %figure; plot(0:h:T,Y); return % TEMP only for testing

    morerror = zeros(4,numel(epsilon));
    redorder = zeros(4,numel(epsilon));
    preerror = zeros(1,numel(epsilon));
    inderror = zeros(1,numel(epsilon));

    ODE = @(f,g,t,x0,u,p) shunt(A,B,C,E,t,x0,u,p,f(0,0,0,0));			% set emgr ODE solver

    rand('seed',1009);
    tic; WC = emgr(F,g,[M,N,Q],[dt,Tf],'c',0,[0,0,0,0,0,0,0,0,0,0,0,0],'r'); toc	% empirical controllability Gramian
    tic; WO = emgr(G,g,[Q,N,M],[dt,Tf],'c',0,[0,0,0,0,0,0,0,0,0,0,0,0],'r'); toc	% empirical observability Gramian
    tic; WX = emgr(F,G,[M,N,Q],[dt,Tf],'y',0,[0,0,0,0,0,0,1,0,0,0,0,0],'r'); toc	% empirical cross Gramian

    obc = ( norm(sum(B,2),2) * norm(sum(C,1),2) );				% Constant error indicator factor

    c = struct('nSnapshots',N,'nModes',0,'tNode',0);

    w = 0.1;									% HAPOD relaxation

    for k = epsilon

        disp(k);
        curr = find(epsilon==k);

        % Cross-Gramian-Based Dominant Subspaces
        [ux,dx1] = hapod_wrapper(WX,k,w,0);
        [vx,dx2] = hapod_wrapper(WX,k,w,1);
        [U,~] = hapod({ux.*dx1,vx.*dx2},k*sqrt(N),'dist_r',w,{c,c});
        V = U;

        y = rk1i(V'*A*U,V'*B,C*U,V'*E*U,[dt,Tf],zeros(size(U,2),1),ut,0,0);	% reduced solution
        morerror(4,curr) = norm(Y(:) - y(:),2)/n2;				% model reduction error
        redorder(4,curr) = size(U,2);						% reduced order
        preerror(curr) = sqrt(obc * k)/n2;					% predicted error
        inderror(curr) = sqrt(obc * norm(WX - ux*(ux'*WX),'fro')/sqrt(N))/n2;	% error indicator

        % Dominant Subspaces Projection Model Reduction
        [uc,dc] = hapod_wrapper(WC,k,w,0);
        [uo,dq] = hapod_wrapper(WO,k,w,0);
        [U,~] = hapod({uc.*sqrt(dc),uo.*sqrt(dq)},k,'dist',w);
        V = U;

        y = rk1i(V'*A*U,V'*B,C*U,V'*E*U,[dt,Tf],zeros(size(U,2),1),ut,0,0);	% reduced solution
        morerror(2,curr) = norm(Y(:) - y(:),2)/n2;				% model reduction error
        redorder(2,curr) = size(U,2);						% reduced order

        % Refined Dominant Subspaces Projection Model Reduction
        fc = sqrt(sum(dc.*dc));
        fo = sqrt(sum(dq.*dq));
        [U,~] = hapod({uc.*(sqrt(dc).*(fo/fc)),uo.*(sqrt(dq))},k,'dist',w);
        V = U;

        y = rk1i(V'*A*U,V'*B,C*U,V'*E*U,[dt,Tf],zeros(size(U,2),1),ut,0,0);	% reduced solution
        morerror(3,curr) = norm(Y(:) - y(:),2)/n2;				% model reduction error
        redorder(3,curr) = size(U,2);						% reduced order

        % Balanced Truncation
        lc = uc.*sqrt(dc);
        lo = uo.*sqrt(dq);
        [ub,db,vb] = svd(lc'*lo,'econ');
        hs = 1.0./sqrt(diag(db))';
        U = lc * ub .* hs;
        V = lo * vb .* hs;
        y = rk1i(V'*A*U,V'*B,C*U,V'*E*U,[dt,Tf],zeros(size(U,2),1),ut,0,0);	% reduced solution
        morerror(1,curr) = norm(Y(:) - y(:),2)/n2;				% model reduction error
        redorder(1,curr) = size(U,2);						% reduced order
    end

    % Plot: Projection Error vs Model Reduction Error
    figure;
    loglog(epsilon,morerror(1,:)','*-','LineWidth',3,'MarkerSize',10);
    hold on;
    loglog(epsilon,morerror(2,:)','*-','LineWidth',3,'MarkerSize',10);
    loglog(epsilon,morerror(3,:)','*--','LineWidth',3,'MarkerSize',10);
    loglog(epsilon,morerror(4,:)','*-','LineWidth',3,'MarkerSize',10);
    loglog(epsilon,preerror,'*-','LineWidth',3,'MarkerSize',10);
    loglog(epsilon,inderror,'*-','LineWidth',3,'MarkerSize',10);
    if(experiment>1)
        loglog(epsilon(3:4),[1;0.1],'k','LineWidth',3);
        loglog(epsilon(3:4),1*ones(2,1),'k','LineWidth',3);
        loglog(ones(2,1)*epsilon(4),[1;0.1],'k','LineWidth',3);
        set(gca,'Xtick',[1e-8,1e-6,1e-4,1e-2]);
    else
        loglog(epsilon(5:6),[1;1e-1],'k','LineWidth',3);
        loglog(epsilon(5:6),ones(2,1),'k','LineWidth',3);
        loglog(ones(2,1)*epsilon(6),[1;1e-1],'k','LineWidth',3);
        set(gca,'Xtick',[1e-11,1e-9,1e-7,1e-5,1e-3]);
    end
    hold off;
    if(experiment>1), xlim([1e-9,1e-1]); ylim([1e-12,1e1]); else, xlim([1e-13,1e-2]); ylim([1e-15,1e1]); end;
    set(gca,'Xdir','reverse','XGrid','on');
    l = legend('LREBT','DSPMR','DSPMR-R','WXDS','H2PRE','H2IND','location','southwest');
    xlabel('Projection Error');
    ylabel('Relative L_2 Model Reduction Error');
    colormap(lines);
    set([gca; findall(gca,'Type','text')],'FontSize',12);
    set(l,'FontSize',12);

    % Plot: Reduced Order vs Model Reduction Error
    figure;
    semilogy(redorder(1,:),morerror(1,:)','*-','LineWidth',3,'MarkerSize',10);
    hold on;
    semilogy(redorder(2,:),morerror(2,:)','*-','LineWidth',3,'MarkerSize',10);
    semilogy(redorder(3,:),morerror(3,:)','*--','LineWidth',3,'MarkerSize',10);
    semilogy(redorder(4,:),morerror(4,:)','*-','LineWidth',3,'MarkerSize',10);
    semilogy(redorder(4,:),preerror,'*-','LineWidth',3,'MarkerSize',10);
    semilogy(redorder(4,:),inderror,'*-','LineWidth',3,'MarkerSize',10);
    set(gca,'XGrid','on');
    hold off;
    if(experiment>1), xlim([0,30]); ylim([1e-12,1e1]); else, xlim([0,35]); ylim([1e-15,1e1]); end;
    l = legend('LREBT','DSPMR','DSPMR-R','WXDS','H2PRE','H2IND','location','northeast');
    xlabel('Reduced Order');
    ylabel('Relative L_2 Model Reduction Error');
    colormap(lines);
    set([gca; findall(gca,'Type','text')],'FontSize',12);
    set(l,'FontSize',12);

    save(['numex',num2str(experiment),'.mat'],'epsilon','morerror','redorder','preerror','inderror');
end

function y = shunt(A,B,C,E,t,x0,u,p,s)

    if(isequal(s,0))
        y = rk1i(A,B,1,E,t,x0,u,p,0);
    else
        y = rk1i(A,C,1,E,t,x0,u,p,1);
    end
end

function y = rk1i(A,B,C,E,t,x0,u,p,s)

    n = numel(x0);
    h = t(1);
    K = floor(t(2)/h) + 1;

    if(s)
        ea = speye(n)-h*(A/E)';
        eb = (B/E)';
    else
        ea = speye(n)-h*(E\A);
        eb = E\B;
    end

    [lea,uea,pea] = lu(ea,'vector');
    eb = eb(pea,:);

    x = x0;
    y = C*x;
    y(:,K) = 0;

    for k = 2:K
        t = (k-1.5)*h;
        x = uea\(lea\(x(pea) + h*eb*u(t) ) );
        y(:,k) = C*x;
    end
end

function [u1,d1] = hapod_wrapper(S,pe,w,rc)

    psize = 200;
    n = ceil(size(S,2)/psize);

    if(rc)
        chunk = @(l) S((l-1)*psize+1:min(size(S,1),l*psize),:)';
    else
        chunk = @(l) S(:,(l-1)*psize+1:min(size(S,2),l*psize));
    end

    u1 = [];
    c1.nLevels = n;

    for k=1:n-1

        [u1,d1,c1] = hapod({chunk(k),u1},pe,'incr_1',w,c1);
    end

    [u1,d1,c1] = hapod({chunk(n),u1},pe,'incr_r',w,c1);
end
