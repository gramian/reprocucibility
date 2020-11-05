function RUNME()
% RUNME (test inverse sylvester procedure)
% by Christian Himpe, 2016 ( http://gramian.de )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%*
    if(exist('emgr')~=2)
        error('emgr not found! Get emgr at: http://gramian.de');
    else
        global ODE; ODE = [];
        fprintf('emgr (version: %g)\n',emgr('version'));
    end

%% SETUP
    J = 1;
    N = 1000;
    O = J;
    [A,B,C] = isp(N,J,1009); A = 0.1*A;
    T = [0.01,1.0];
    L = floor(T(2)/T(1)) + 1;
    U = [ones(J,1)./T(1),zeros(J,L-1)];
    randn('seed',1009);
    V = randn(J,L);
    X = zeros(N,1);

    LIN = @(x,u,p) A*x + B*u;
    OUT = @(x,u,p) C*x;

    %figure; plot(0:T(1):T(2),ODE(LIN,OUT,T,X,U,0)); return;

%% MAIN

    m = 100; % N

    l1 = zeros(m,3);
    l2 = zeros(m,3);
    l3 = zeros(m,3);
    h2 = zeros(m,3);
    h8 = zeros(m,3);

    tic;
    WX1 = sylvester(A,A,-B*C);
    [UU,D,VV] = svd(WX1); D = diag(D); VV = VV';
    OFFLINE_WX1 = toc
    [l1(:,1),l2(:,1),l8(:,1),h2(:,1),h8(:,1)] = assess(ODE,LIN,OUT,T,X,V,0,UU,D,VV,VV*B,C*UU);

    tic;
    ADJ = @(x,u,p) A'*x + C'*u;
    WX2 = emgr(LIN,ADJ,[J,N,O],T,'y');
    [UU,D,VV] = svd(WX2); D = diag(D); VV = VV';
    OFFLINE_WX2 = toc
    [l1(:,2),l2(:,2),l8(:,2),h2(:,2),h8(:,2)] = assess(ODE,LIN,OUT,T,X,V,0,UU,D,VV,VV*B,C*UU);

    tic;
    WX3 = emgr(LIN,OUT,[J,N,O],T,'x');
    [UU,D,VV] = svd(WX3); D = diag(D); VV = UU';
    OFFLINE_WX3 = toc
    [l1(:,3),l2(:,3),l8(:,3),h2(:,3),h8(:,3)] = assess(ODE,LIN,OUT,T,X,V,0,UU,D,VV,VV*B,C*UU);

    visual(l1,'l1');
    visual(l2,'l2');
    visual(l8,'l8');
    visual(h2,'h2');
    visual(h8,'h8');
end

%% ======== Visual ========
function visual(d,n)

    m = 100; %size(d,1);

    figure('Name',n,'NumberTitle','off');
    semilogy(1:m,d(:,1),'r','linewidth',10); hold on;
    semilogy(1:m,d(:,2),'g','linewidth',10);
    semilogy(1:m,d(:,3),'b--','linewidth',10); hold off;
    xlim([1,m]);
    ylim([1e-18,1]);
    pbaspect([1,0.5,1]);
    set([gca; findall(gca,'Type','text')],'FontSize',14);
    l = legend('WX1','WX2','WX3','location','northeast');
    set(l,'FontSize',10);
    print('-depsc',[n,'.eps']);
end

%% ======== Assess ========
function [l1,l2,l8,h2,h8] = assess(ODE,LIN,OUT,T,X,U,P,UU,D,VV,BB,CC)

    M = 100; %numel(X);

    Y = ODE(LIN,OUT,T,X,U,P);
    L1 = norm(Y(:),1);
    L2 = norm(Y(:),2);
    L8 = norm(Y(:),Inf);
    H2 = abs(sqrt(trace(CC*diag(D)*BB)));
    H8 = 2.0*sum(D);

    l1 = zeros(M,1);
    l2 = zeros(M,1);
    l8 = zeros(M,1);
    h2 = zeros(M,1);
    h8 = zeros(M,1);

    for I=1:M
        uu = UU(:,1:I);
        vv = VV(1:I,:);
        lin = @(x,u,p) vv*LIN(uu*x,u,p);
        out = @(x,u,p) OUT(uu*x,u,p);
        y = ODE(lin,out,T,vv*X,U,P);
        l1(I) = norm(Y(:)-y(:),1)/L1;
        l2(I) = norm(Y(:)-y(:),2)/L2;
        l8(I) = norm(Y(:)-y(:),Inf)/L8;
        h2(I) = abs(sqrt(trace(CC(:,I+1:end)*diag(D(I+1:end))*BB(I+1:end,:))))/H2;
        h8(I) = 2.0*sum(D(I+1:end))/H8;
    end;
end
