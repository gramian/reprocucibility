function nonsym(e)
% nonsym (Non-Symmetric Cross Gramian Companion Code)
% by Christian Himpe, 2014-2016 ( http://gramian.de )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%*
    if(exist('emgr')~=2)
        disp('emgr framework is required. Download at http://gramian.de/emgr.m');
        return;
    else
        global ODE;
        fprintf('emgr (version: %g)\n',emgr('version'));
    end

    T = [0.002,1.0];
    L = floor(T(2))/T(1) + 1;

    switch(e)

        case 0,
            J = 8;
            N = 1024;
            O = J;
            U = [ones(J,1),zeros(J,L-1)];
            X = zeros(N,1);            

            rand('seed',1009);
            A = -0.02*gallery('lehmer',N);
            B = rand(N,J);
            C = B';

            LIN = @(x,u,p) A*x + B*u;
            OUT = @(x,u,p) C*x;

            b = sum(B,2);
            c = sum(C,1);

            TLIN = @(x,u,p) A*x + b*u;
            TOUT = @(x,u,p) c*x;

            tic;
            WX = emgr(TLIN,TOUT,[1,N,1],T,'x');
            OFFLINE_WX0 = toc

            tic;
            WZ = emgr(LIN,OUT,[J,N,O],T,'x',0,[0,0,0,0,0,0,1,0,0,0,0,0]);
            OFFLINE_WZ0 = toc

            FROBNORM = norm(WX-WZ,'fro')

        case 1,
            J = 8;
            N = 1024;
            O = J;
            U = [ones(J,1),zeros(J,L-1)];            
            X = zeros(N,1);

            rand('seed',1009);
            A = -0.02*gallery('lehmer',N);
            B = rand(N,J);
            C = B';

            LIN = @(x,u,p) A*x + B*u;
            OUT = @(x,u,p) C*x;

            Y = ODE(LIN,OUT,T,X,U,0);

            tic;
            WC = emgr(LIN,OUT,[J,N,O],T,'c');
            WO = emgr(LIN,OUT,[J,N,O],T,'o');
            [VV,D,UU] = balance(WC,WO);
            OFFLINE_BT1 = toc
            l2bt1 = errortest(A,B,C,X,U,T,UU,VV,Y);

            tic;
            WX = emgr(LIN,OUT,[J,N,O],T,'x');
            [UU,D,VV] = svd(WX);
            OFFLINE_WX1 = toc
            l2wx1 = errortest(A,B,C,X,U,T,UU,UU',Y);

            tic;
            WZ = emgr(LIN,OUT,[J,N,O],T,'x',0,[0,0,0,0,0,0,1,0,0,0,0,0]);
            [UU,D,VV] = svd(WZ);
            OFFLINE_WZ1 = toc
            l2wz1 = errortest(A,B,C,X,U,T,UU,UU',Y);

            makeplot({l2wz1,l2wx1,l2bt1},{'Nonsym. Cross Gramian','Cross Gramian','Balanced Truncation'},'figure1');

        case 2,
            J = 1;
            N = 1024;
            O = 8;
            U = [ones(J,1),zeros(J,L-1)];            
            X = zeros(N,1);

            rand('seed',1009);
            A = -0.02*gallery('lehmer',N);
            k = 2*rand(N,1);
            K = diag(k);
            KI = diag(1.0/k);
            A = A*K;
            B = rand(N,J);
            C = rand(O,N);

            LIN = @(x,u,p) A*x + B*u;
            OUT = @(x,u,p) C*x;

            Y = ODE(LIN,OUT,T,X,U,0);

            LINA = @(x,u,p) A*x + [K*C',B]*u;
            OUTA = @(x,u,p) [C;B'*KI]*x;

            tic;
            WC = emgr(LIN,OUT,[J,N,O],T,'c');
            WO = emgr(LIN,OUT,[J,N,O],T,'o');
            [VV,D,UU] = balance(WC,WO);
            OFFLINE_BT2 = toc
            l2bt2 = errortest(A,B,C,X,U,T,UU,VV,Y);

            tic;
            WX = emgr(LINA,OUTA,[J+O,N,J+O],T,'x');
            [UU,D,VV] = svd(WX);
            OFFLINE_WX2 = toc
            l2wx2 = errortest(A,B,C,X,U,T,UU,UU',Y);

            tic;
            WZ = emgr(LIN,OUT,[J,N,O],T,'x',0,[0,0,0,0,0,0,1,0,0,0,0,0]);
            [UU,D,VV] = svd(WZ);
            OFFLINE_WZ2 = toc
            l2wz2 = errortest(A,B,C,X,U,T,UU,UU',Y);

            makeplot({l2wz2,l2wx2,l2bt2},{'Nonsym. Cross Gramian','Approx. Cross Gramian','Balanced Truncation'},'figure2');

        case 3,
            J = 8;
            N = 1024;
            O = 8;
            U = [ones(J,1),zeros(J,L-1)];         
            X = zeros(N,1);

            rand('seed',1009);
            A = -0.02*gallery('lehmer',N);
            B = rand(N,J);
            C = rand(O,N);

            LIN = @(x,u,p) A*x + B*u;
            OUT = @(x,u,p) C*x;

            Y = ODE(LIN,OUT,T,X,U,0);

            LINA = @(x,u,p) A*x + B*u;
            OUTA = @(x,u,p) C*x;

            tic;
            WC = emgr(LIN,OUT,[J,N,O],T,'c');
            WO = emgr(LIN,OUT,[J,N,O],T,'o');
            [VV,D,UU] = balance(WC,WO);
            OFFLINE_BT3 = toc
            l2bt3 = errortest(A,B,C,X,U,T,UU,VV,Y);

            tic;
            WX = emgr(LIN,OUT,[J,N,O],T,'x');
            [UU,D,VV] = svd(WX);
            OFFLINE_WX3 = toc
            l2wx3 = errortest(A,B,C,X,U,T,UU,UU',Y);

            tic;
            WZ = emgr(LIN,OUT,[J,N,O],T,'x',0,[0,0,0,0,0,0,1,0,0,0,0,0]);
            [UU,D,VV] = svd(WZ);
            OFFLINE_WZ3 = toc
            l2wz3 = errortest(A,B,C,X,U,T,UU,UU',Y);

            makeplot({l2wz3,l2wx3,l2bt3},{'Nonsym. Cross Gramian','Cross Gramian','Balanced Truncation'},'figure3');

        case 4,
            if(exist('/tmp/iss.mat','file')==0)
                urlwrite('http://slicot.org/objects/software/shared/bench-data/iss.zip','/tmp/iss.zip');
                unzip('/tmp/iss.zip','/tmp');
            end
            load(['/tmp/iss.mat']);
            N = size(A,2);
            n = N/2;
            J = size(B,2);
            O = size(C,1);
            T = [0.01,1.0];
            L = floor(T(2)/T(1)) + 1;
            U = [ones(J,1),zeros(J,L-1)];
            X = zeros(N,1);

            LIN = @(x,u,p) A*x + B*u;
            OUT = @(x,u,p) C*x;

            Y = ODE(LIN,OUT,T,X,U,0);
            %figure; plot(0:T(1):T(2),Y); return;

            tic;
            WZ = emgr(LIN,OUT,[J,N,O],T,'x',0,[0,0,0,0,0,0,1,0,0,0,0,0]);
            [UU,D,VV] = svd(WZ(n+1:N,n+1:N));
            OFFLINE_WZ4 = toc
            l2wz4 = errortest_so(A,B,C,X,U,T,UU,UU',Y);

            tic;
            WX = emgr(LIN,OUT,[J,N,O],T,'x',0);
            [UU,D,VV] = svd(WX(n+1:N,n+1:N));
            OFFLINE_WX4 = toc
            l2wx4 = errortest_so(A,B,C,X,U,T,UU,UU',Y);

            tic;
            WC = emgr(LIN,OUT,[J,N,O],T,'c');
            WO = emgr(LIN,OUT,[J,N,O],T,'o');
            [UU,D,VV] = balance(WC(n+1:N,n+1:N),WO(n+1:N,n+1:N));
            OFFLINE_BT4 = toc
            l2bt4 = errortest_so(A,B,C,X,U,T,VV,UU,Y);

            makeplot({l2wz4,l2wx4,l2bt4},{'Nonsym Cross Gramian','Cross Gramian','Balanced Truncation'},'figure4');

        case 5,
            J = 1;
            N = 1006;
            O = 8;
            T = [0.001,0.1];
            L = floor(T(2)/T(1)) + 1;
            U = [ones(J,1),zeros(J,L-1)];         
            X = zeros(N,1);

            rand('seed',1009);
            A = zeros(N,N);
            A(1:2,1:2) = [-1,100;-100,-1];
            A(3:4,3:4) = [-1,200;-200,-1];
            A(5:6,5:6) = [-1,400;-400,-1];
            A(7:end,7:end) = -diag(linspace(1,N-6,N-6));
            B = [10*ones(6,1);ones(N-6,1)];
            C = [10*rand(6,O);rand(N-6,O)]';

            LIN = @(x,u,p) A*x + B*u;
            OUT = @(x,u,p) C*x;

            Y = ODE(LIN,OUT,T,X,U,0);
            %figure; plot(0:T(1):T(2),Y); return;

            tic;
            WZ = emgr(LIN,OUT,[J,N,O],T,'x',0,[0,0,0,0,0,0,1,0,0,0,0,0]);
            [UU,D,VV] = svd(WZ);
            OFFLINE_WZ5 = toc
            l2wz5 = errortest(A,B,C,X,U,T,UU,UU',Y);

            tic;
            WC = emgr(LIN,OUT,[J,N,O],T,'c');
            WO = emgr(LIN,OUT,[J,N,O],T,'o');
            [UU,D,VV] = balance(WC,WO);
            OFFLINE_BT5 = toc
            l2bt5 = errortest(A,B,C,X,U,T,VV,UU,Y);

            makeplot({l2wz5,l2bt5},{'Nonsym Cross Gramian','Balanced Truncation'},'figure5');
    end
end

%% ======== Errortest ========
function e = errortest(A,B,C,X,U,T,UU,VV,Y)

    global ODE;

    N = size(A,1);

    e = zeros(1,N-1);

    for I=1:N-1
        uu = UU(:,1:I);
        vv = VV(1:I,:);
        a = vv*A*uu;
        b = vv*B;
        c = C*uu;
        x = vv*X;
        lin = @(x,u,p) a*x + b*u;
        out = @(x,u,p) c*x;
        y = ODE(lin,out,T,x,U,0); % Reduced Order
        e(I) = norm(Y(:)-y(:),2)/norm(Y(:),2);
    end
end

%% ======== 2nd Order Errortest ========
function e = errortest_so(A,B,C,X,U,T,UU,VV,Y)

    global ODE;

    N = size(A,1);
    n = N/2;
    J = size(B,2);
    O = size(C,1);

    e = zeros(1,n-1);

    for I=1:n-1
        uu = UU(:,1:I);
        vv = VV(1:I,:);
        a = [zeros(I,I),eye(I);vv*A(n+1:N,1:n)*uu,vv*A(n+1:N,n+1:N)*uu];
        b = [zeros(I,J);vv*B(n+1:N,:)];
        c = [zeros(O,I),C(:,n+1:N)*uu];
        x = zeros(2*I,1);
        lin = @(x,u,p) a*x + b*u;
        out = @(x,u,p) c*x;
        y = ODE(lin,out,T,x,U,0); % Reduced Order
        e(I) = norm(Y(:)-y(:),2)/norm(Y(:),2);
    end
end

%% ======== Make Plots ========
function makeplot(e,l,f)

    n = numel(e{1});
    m = 100;

    figure;
    semilogy(1:m,e{1}(1:m),'r','linewidth',4);
    hold on;
    semilogy(1:m,e{2}(1:m),'g','linewidth',4);
    if(numel(e)>=3), semilogy(1:m,e{3}(1:m),'b','linewidth',4); end;
    if(numel(e)>=4), semilogy(1:m,e{4}(1:m),'k','linewidth',4); end;
    hold off;

    xlim([1,m]);
    ylim([1e-16,1]);
    pbaspect([2,1,1]);
    xlabel('Reduced State Space Dimension');
    ylabel('Relative Error');
    if(numel(e)==4)
        legend(l{1},l{2},l{3},l{4},'location','northeast');
    elseif(numel(e)==3)
        legend(l{1},l{2},l{3},'location','northeast');
    else
        legend(l{1},l{2},'location','northeast');
    end

    print('-dpdf',[f,'.pdf']);
end


%% ======== Balancer ========
function [X,Y,Z] = balance(WC,WO)

    [L,D,l] = svd(WC); LC = L*diag(sqrt(diag(D)));
    [L,D,l] = svd(WO); LO = L*diag(sqrt(diag(D)));
    [U,Y,V] = svd(LO'*LC);
    X = ( LO*U*diag(1.0./sqrt(diag(Y))) )';
    Z =   LC*V*diag(1.0./sqrt(diag(Y)));
end

