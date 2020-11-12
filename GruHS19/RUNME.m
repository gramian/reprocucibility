function RUNME()
%%% summary: linear advection example
%%% project: GAMM / PAMM 2019
%%% authors: Christian Himpe ( 0000-0003-2194-6754 )
%%% license: BSD-2-Clause (2019)

    global ODE;
    global STAGES;
    STAGES = 1;
    EMGR_VERSION = emgr('version')

    % System Dimensions
    M = 1;
    N = 1000;
    Q = 1;

    % Time Discretization
    dt = 0.5/N;
    Tf = 1.5;

    % System Matrices
    A = spdiags(N*ones(N,1)*[1,-1],[-1,0],N,N);
    B = sparse(1,1,N,N,1);
    C = sparse(1,N,1.0,1,N);

    % Zero Initial Condition
    X0 = zeros(N,1);

    % Test Input
    vt = @(t) exp(-100*(t-0.2).^2);
    %figure; plot([0:dt:Tf],vt([0:dt:Tf])); return; % ONLY for visualization of test input

    % System Setup
    f = @(x,u,p,t) A*x + B*u;
    g = @(x,u,p,t) C*x;
    F = @(x,u,p,t) A'*x + C'*u;

    %X = ODE(f,@(x,u,p,t) x,[dt,Tf],X0,vt,0); figure; plot([0:dt:Tf],X); return; % ONLY for visualization of test state evolution

    if not(exist('OCTAVE_VERSION','builtin')) % MATLAB requires dense matrices for sylvester solver

        A = full(A);
        B = full(B);
        C = full(C);
    end

    % Compute Algebraic System Gramians
    WCl = sylvester(A,A',-B*B'); % Controllability Gramian
    WOl = sylvester(A',A,-C'*C); % Observability Gramian
    WYl = sylvester(A,A,-B*C);   % Cross Gramian

    % Compute Empirical System Gramians
    WCn = emgr(f,g,[M,N,Q],[dt,Tf],'c',0,[0,0,0,0,0,0,0,0,0,0,0,0],"s"); % Empirical Controllability Gramian
    WOn = emgr(F,g,[M,N,Q],[dt,Tf],'c',0,[0,0,0,0,0,0,0,0,0,0,0,0],"s"); % Empirical (Linear) Observability Gramian
    WYn = emgr(f,F,[M,N,Q],[dt,Tf],'y',0,[0,0,0,0,0,0,0,0,0,0,0,0],"s"); % Empirical (Linear) Cross Gramian

    % Proper Orthogonal Decomposition (POD)
    [U_PODl,D_PODl,~] = svd(WCl);
    [U_PODn,D_PODn,~] = svd(WCn);

    % Balanced Truncation (BT)
    [V_BTl,D_BTl,U_BTl] = balco(WCl,WOl);
    [V_BTn,D_BTn,U_BTn] = balco(WCn,WOn);

    % Dominant Subspaces Projection Model Reduction (DSPMR)
    [U_DSl,D_DSl,V_DSl] = svd(WYl); U_DSl = U_DSl*D_DSl; V_DSl = V_DSl*D_DSl;
    [U_DSn,D_DSn,V_DSn] = svd(WYn); U_DSn = U_DSn*D_DSn; V_DSn = V_DSn*D_DSn;

    % Algebraic Reductors
    red{1} = @(k) deal(U_PODl(:,1:k),U_PODl(:,1:k));
    red{2} = @(k) deal(domsub(U_DSl(:,1:k),V_DSl(:,1:k),k));
    red{3} = @(k) deal(U_BTl(:,1:k),V_BTl(:,1:k));

    % Empirical Reductors
    red{4} = @(k) deal(U_PODn(:,1:k),U_PODn(:,1:k));
    red{5} = @(k) deal(domsub(U_DSn(:,1:k),V_DSn(:,1:k),k));
    red{6} = @(k) deal(U_BTn(:,1:k),V_BTn(:,1:k));

    % Test Solution
    y = ODE(f,g,[dt,Tf],X0,vt,0);
    ny2 = norm(y,2);

    % Test Reduced Orders
    r = round(linspace(1,500,100));

    % Reduced Order Models
    for k = r
        for l = [1,2,3,4,5,6]

            [ur,vr] = red{l}(k);
            a = vr'*A*ur;
            b = vr'*B;
            c = C*ur;
            fr = @(x,u,p,t) a*x + b*u;
            gr = @(x,u,p,t) c*x;
            yr = ODE(fr,gr,[dt,Tf],vr'*X0,vt,0);
            l2(l,find(r==k)) = norm(y-yr,2)/ny2;
        end
    end

    % Visualization
    figure;
    semilogy(r,l2(1,:),'Linewidth',3);
    hold on;
    semilogy(r,l2(2,:),'Linewidth',3);
    semilogy(r,l2(3,:),'Linewidth',3);
    semilogy(r,l2(4,:),'Linewidth',3);
    semilogy(r,l2(5,:),'Linewidth',3);
    semilogy(r,l2(6,:),'Linewidth',3);
    hold off;
    xlim([1,500]);
    ylim([1e-15,1]);
    set(gca,'XTick',[1,50:50:500]);
    set(gca,'YTick',[1e-15,1e-12,1e-9,1e-6,1e-3,1e-0]);
    xlabel('State Dimension');
    ylabel('Relative Output Error');
    colormap(lines);
    pbaspect([2,1,1]);
    legend('POD','BT','DSPMR','Empirical POD','Empirical BT','Empirical DSPMR');
end

function [S,T] = domsub(U,V,k)
% Dominant Subspace

    [t,~,~] = svd([U,V],'econ');
    S = t(:,1:k);
    T = S;
end

function [TR,HSV,TL] = balco(WC,WO)
% Balancing Transformation

    [LC,EC] = eig(WC,'vector');
    [ec,IN] = sort(sqrt(abs(EC)),'descend');
    [LO,EO] = eig(WO,'vector');
    [eo,IN] = sort(sqrt(abs(EO)),'descend');
    LC = LC(:,IN)*diag(ec);
    LO = LO(:,IN)*diag(eo);

    [U,HSV,V] = svd(LC'*LO);
    HSV = diag(HSV);
    y = sqrt(1.0./(HSV + 1e-14));
    TR = LO*V*diag(y);
    TL = LC*U*diag(y);
end
