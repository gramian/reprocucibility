function mpe(n)
% mpie (Companion Code)
% by Christian Himpe, 2014 ( http://gramian.de )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%*

if(exist('emgr')~=2) disp('emgr framework is required. Download at http://gramian.de/emgr.m'); return; end

rand('seed',1009);
randn('seed',1009);

%%%%%%%% Setup %%%%%%%%

if(nargin<1), n = 16; end;

J = n;
N = J*J;
O = J;
M = 100;

S = 0;
h = 0.01;
T = 1.0;

A = sprand(N,N,0.1)-N*diag(0.25+0.5*rand(N,1));
A = 0.5*(A+A');
B = rand(N,J);
C = B';
X = zeros(N,1);
U = [ones(J,1),zeros(J,(T/h)-1)];
P = 0.1*rand(N,1);

IS_STABLE = all(eig(A)<0)

Fl = @(x,u,p) A*x + B*u + p;
Fa = @(x,u,p) A'*x + C'*u;
Fn = @(x,u,p) A*tanh(0.25*x) + B*u + p;
Gl = @(x,u,p) C*x;

gn = @(x) exp(x)+x-1.0;
A1 = spdiags(ones(M-1,1),-1,M,M)-speye(M);
A2 = spdiags([ones(M-1,1);0],0,M,M)-spdiags(ones(M,1),1,M,M);
Hn = @(x,u,p) gn(A1*x)-gn(A2*x) + [u;sparse(M-1,1)];
Gn = @(x,u,p) x(1);
Xn = zeros(M,1);
Un = exp(-[h:h:T]);

l2norm = @(m) sqrt(sum(sum(m.*m,1)));
l8norm = @(m) max(sqrt(sum(m.*m,1)));

OFFLINE = zeros(1,20);
SINGVAL = zeros(20,N);
L2ERROR = zeros(20,N-1);
L8ERROR = zeros(20,N-1);
l2error = zeros(N-1,N-1,20);
l8error = zeros(N-1,N-1,20);

Yl = assess(Fl,Gl,X,U,P,S,h,T,speye(N),speye(N),N,speye(N),speye(N),N);
Yn = assess(Fn,Gl,X,U,P,S,h,T,speye(N),speye(N),N,speye(N),speye(N),N);
Bn = assess(Hn,Gn,Xn,Un,0,S,h,T,speye(M),speye(M),M,speye(M),speye(M),M);

IND = [1,J:J:N];

%subplot(1,2,1); plot(h:h:T,Yl'); subplot(1,2,2); plot(h:h:T,Yn'); return;

%%%%%%%% Linear State Reduction %%%%%%%%

% BPOD
tic;
WY = emgr(Fl,Fa,[J N O],[S h T],'y',P,0,U);
[UY DY VY] = svd(WY);
OFFLINE(1) = toc;
SINGVAL(1,:) = diag(DY);
for I=1:N-1, y = assess(Fl,Gl,X,U,P,S,h,T,UY',UY,I,speye(N),speye(N),N); L2ERROR(1,I) = l2norm(Yl-y); L8ERROR(1,I) = l8norm(Yl-y); end;
fprintf('#');

% BT
tic;
WC = emgr(Fl,Gl,[J N O],[S h T],'c',P,0,U);
twc = toc;
WO = emgr(Fl,Gl,[J N O],[S h T],'o',P,0,U);
two = toc - twc;
[UB DB VB] = bal(WC,WO);
OFFLINE(2) = toc;
SINGVLA(2,:) = diag(DB);
for I=1:N-1, y = assess(Fl,Gl,X,U,P,S,h,T,UB,VB,I,speye(N),speye(N),N); L2ERROR(2,I) = l2norm(Yl-y); L8ERROR(2,I) = l8norm(Yl-y); end;
fprintf('#');

% WX
tic;
WX = emgr(Fl,Gl,[J N O],[S h T],'x',P,0,U);
[UX DX VX] = svd(WX);
OFFLINE(3) = toc;
SINGVAL(3,:) = diag(DX);
for I=1:N-1, y = assess(Fl,Gl,X,U,P,S,h,T,UX',UX,I,speye(N),speye(N),N); L2ERROR(3,I) = l2norm(Yl-y); L8ERROR(3,I) = l8norm(Yl-y); end;
fprintf('#');

%%%%%%%% Linear Parameter Reduction %%%%%%%%

% WS
tic;
WS = emgr(Fl,Gl,[J N O],[S h T],'s',P,0,U,1);
[US DS VS] = svd(WS{2});
OFFLINE(4) = toc;
SINGVAL(4,:) = diag(DS);
for I=1:N-1, y = assess(Fl,Gl,X,U,P,S,h,T,speye(N),speye(N),N,US',US,I); L2ERROR(4,I) = l2norm(Yl-y); L8ERROR(4,I) = l8norm(Yl-y); end;
fprintf('#');

% WI
tic;
WI = emgr(Fl,Gl,[J N O],[S h T],'i',P,0,U,0,1);
[UI DI VI] = svd(WI{2});
OFFLINE(5) = toc;
SINGVAL(5,:) = diag(DI);
for I=1:N-1, y = assess(Fl,Gl,X,U,P,S,h,T,speye(N),speye(N),N,UI',UI,I); L2ERROR(5,I) = l2norm(Yl-y); L8ERROR(5,I) = l8norm(Yl-y); end;
fprintf('#');

% WJ
tic;
WJ = emgr(Fl,Gl,[J N O],[S h T],'j',P,0,U,0,1);
[UJ DJ VJ] = svd(WJ{2});
OFFLINE(6) = toc;
SINGVAL(6,:) = diag(DJ);
for I=1:N-1, y = assess(Fl,Gl,X,U,P,S,h,T,speye(N),speye(N),N,UJ',UJ,I); L2ERROR(6,I) = l2norm(Yl-y); L8ERROR(6,I) = l8norm(Yl-y); end;
fprintf('#');

%%%%%%%% Linear Combined Reduction %%%%%%%%

% WS+WO
tic;
[UCC DCC VCC] = bal(WS{1},WO);
OFFLINE(7) = toc + OFFLINE(4) + two;
SINGVAL(7,:) = diag(DCC);
for I=1:J, for K=1:J, y = assess(Fl,Gl,X,U,P,S,h,T,UCC,VCC,IND(I),US',US,IND(K)); l2error(I,K,7) = l2norm(Yl-y); l8error(I,K,7) = l8norm(Yl-y); end; end;
for I=1:N-1, y = assess(Fl,Gl,X,U,P,S,h,T,UCC,VCC,I,US',US,I); L2ERROR(7,I) = l2norm(Yl-y); end;
fprintf('#');

% WC+WI
tic;
[UCO DCO VCO] = bal(WC,WI{1});
OFFLINE(8) = toc + OFFLINE(5) + twc;
SINGVAL(8,:) = diag(DCO);
for I=1:J, for K=1:J, y = assess(Fl,Gl,X,U,P,S,h,T,UCO,VCO,IND(I),UI',UI,IND(K)); l2error(I,K,8) = l2norm(Yl-y); l8error(I,K,8) = l8norm(Yl-y); end; end;
for I=1:N-1, y = assess(Fl,Gl,X,U,P,S,h,T,UCO,VCO,I,UI',UI,I); L2ERROR(8,I) = l2norm(Yl-y); end;
fprintf('#');

% WJ
tic;
[UCX DCX VCX] = svd(WJ{1});
OFFLINE(9) = toc + OFFLINE(6);
SINGVAL(9,:) = diag(DCX);
for I=1:J, for K=1:J, y = assess(Fl,Gl,X,U,P,S,h,T,UCX',UCX,IND(I),UJ',UJ,IND(K)); l2error(I,K,9) = l2norm(Yl-y); l8error(I,K,9) = l8norm(Yl-y); end; end;
for I=1:N-1, y = assess(Fl,Gl,X,U,P,S,h,T,UCX',UCX,I,UJ',UJ,I); L2ERROR(9,I) = l2norm(Yl-y); end;
fprintf('# ');

%%%%%%%% Nonlinear State Reduction %%%%%%%%

% BPOD
tic;
WC = emgr(Fn,Gl,[J N O],[S h T],'c',P,0,U);
twc = toc;
WO = emgr(Fn,Gl,[J N O],[S h T],'o',P,0,U);
two = toc - twc;
[UY DY VY] = svd(WC*WO);
OFFLINE(10) = toc;
SINGVAL(10,:) = diag(DY);
for I=1:N-1, y = assess(Fn,Gl,X,U,P,S,h,T,UY',UY,I,speye(N),speye(N),N); L2ERROR(10,I) = l2norm(Yn-y); L8ERROR(10,I) = l8norm(Yn-y); end;
fprintf('#');

% BT
tic;
[UB DB VB] = bal(WC,WO);
OFFLINE(11) = toc + twc + two;
SINGVAL(11,:) = diag(DB);
for I=1:N-1, y = assess(Fn,Gl,X,U,P,S,h,T,UB,VB,I,speye(N),speye(N),N); L2ERROR(11,I) = l2norm(Yn-y); L8ERROR(11,I) = l8norm(Yn-y); end;
fprintf('#');

% WX
tic;
WX = emgr(Fn,Gl,[J N O],[S h T],'x',P,0,U);
[UX DX UX] = svd(WX);
OFFLINE(12) = toc;
SINGVAL(12,:) = diag(DX);
for I=1:N-1, y = assess(Fn,Gl,X,U,P,S,h,T,UX',UX,I,speye(N),speye(N),N); L2ERROR(12,I) = l2norm(Yn-y); L8ERROR(12,I) = l8norm(Yn-y); end;
fprintf('#');

%%%%%%%% Nonlinear Parameter Reduction %%%%%%%%

% WS
tic;
WS = emgr(Fn,Gl,[J N O],[S h T],'s',P,0,U,1);
[US DS VS] = svd(WS{2});
OFFLINE(13) = toc;
SINGVAL(13,:) = diag(DS);
for I=1:N-1, y = assess(Fn,Gl,X,U,P,S,h,T,speye(N),speye(N),N,US',US,I); L2ERROR(13,I) = l2norm(Yn-y); L8ERROR(13,I) = l8norm(Yn-y); end;
fprintf('#');

% WI
tic;
WI = emgr(Fn,Gl,[J N O],[S h T],'i',P,0,U,0,1);
[UI DI VI] = svd(WI{2});
OFFLINE(14) = toc;
SINGVAL(14,:) = diag(DI);
for I=1:N-1, y = assess(Fn,Gl,X,U,P,S,h,T,speye(N),speye(N),N,UI',UI,I); L2ERROR(14,I) = l2norm(Yn-y); L8ERROR(14,I) = l8norm(Yn-y); end;
fprintf('#');

% WJ
tic;
WJ = emgr(Fn,Gl,[J N O],[S h T],'j',P,0,U,0,1);
[UJ DJ VJ] = svd(WJ{2});
OFFLINE(15) = toc;
SINGVAL(15,:) = diag(DJ);
for I=1:N-1, y = assess(Fn,Gl,X,U,P,S,h,T,speye(N),speye(N),N,UJ',UJ,I); L2ERROR(15,I) = l2norm(Yn-y); L8ERROR(15,I) = l8norm(Yn-y); end;
fprintf('#');

%%%%%%%% Nonlinear Combined Reduction %%%%%%%%

% WS+WO
tic;
[UCC DCC VCC] = bal(WS{1},WO);
OFFLINE(16) = toc + OFFLINE(13) + two;
SINGVAL(16,:) = diag(DCC);
for I=1:J, for K=1:J, y = assess(Fn,Gl,X,U,P,S,h,T,UCC,VCC,IND(I),US',US,IND(K)); l2error(I,K,16) = l2norm(Yn-y); l8error(I,K,16) = l8norm(Yn-y); end; end;
for I=1:N-1, y = assess(Fn,Gl,X,U,P,S,h,T,UCC,VCC,I,US',US,I); L2ERROR(16,I) = l2norm(Yn-y); end;
fprintf('#');

% WC+WI
tic;
[UCO DCO VCO] = bal(WC,WI{1});
OFFLINE(17) = toc + OFFLINE(14) + twc;
SINGVAL(17,:) = diag(DCO);
for I=1:J, for K=1:J, y = assess(Fn,Gl,X,U,P,S,h,T,UCO,VCO,IND(I),UI',UI,IND(K)); l2error(I,K,17) = l2norm(Yn-y); l8error(I,K,17) = l8norm(Yn-y); end; end;
for I=1:N-1, y = assess(Fn,Gl,X,U,P,S,h,T,UCO,VCO,I,UI',UI,I); L2ERROR(17,I) = l2norm(Yn-y); end;
fprintf('#');

% WJ
tic;
[UCX DCX VCX] = svd(WJ{1});
OFFLINE(18) = toc + OFFLINE(15);
SINGVAL(18,:) = diag(DCX);
for I=1:J, for K=1:J, y = assess(Fn,Gl,X,U,P,S,h,T,UCX',UCX,IND(I),UJ',UJ,IND(K)); l2error(I,K,18) = l2norm(Yn-y); l8error(I,K,18) = l8norm(Yn-y); end; end;
for I=1:N-1, y = assess(Fn,Gl,X,U,P,S,h,T,UCX',UCX,I,UJ',UJ,I); L2ERROR(18,I) = l2norm(Yn-y); end;
fprintf('# ');

%%%%%%%% Nonlinear Benchmark State Reduction %%%%%%%%

% BT
tic;
WC = emgr(Hn,Gn,[1 M 1],[S h T],'c',P,0,Un);
WO = emgr(Hn,Gn,[1 M 1],[S h T],'o',P,0,Un);
[UB DB VB] = bal(WC,WO);
OFFLINE(19) = toc;
%SINGVLA(19,:) = diag(DB);
for I=1:M-1, y = assess(Hn,Gn,Xn,Un,0,S,h,T,UB,VB,I,speye(M),speye(M),M); L2ERROR(19,I) = l2norm(Bn-y); L8ERROR(19,I) = l8norm(Bn-y); end;
fprintf('#');

% WX
tic;
WX = emgr(Hn,Gn,[1 M 1],[S h T],'x',P,0,Un);
[UX DX VX] = svd(WX);
OFFLINE(20) = toc;
%SINGVAL(20,:) = diag(DX);
for I=1:M-1, y = assess(Hn,Gn,Xn,Un,0,S,h,T,UX',UX,I,speye(M),speye(M),M); L2ERROR(20,I) = l2norm(Bn-y); L8ERROR(20,I) = l8norm(Bn-y); end;
fprintf('#\n');

%%%% %%%% %%%% %%%%

ORIGLL2 = l2norm(Yl);
ORIGLL8 = l8norm(Yl);
ORIGNL2 = l2norm(Yn);
ORIGNL8 = l8norm(Yn);
ORIGBL2 = l2norm(Bn);
ORIGBL8 = l8norm(Bn);

save('mpe.mat','N','M','OFFLINE','SINGVAL','L2ERROR','L8ERROR','l2error','l8error','ORIGLL2','ORIGLL8','ORIGNL2','ORIGNL8','ORIGBL2','ORIGBL8');


end

%%%%

function [U Y V] = bal(WC,WO)

    [L D l] = svd(WC); LC = L*diag(sqrt(diag(D)));
    [L D l] = svd(WO); LO = L*diag(sqrt(diag(D)));
    [X Y Z] = svd(LO'*LC);
    U = ( LO*X*diag(1.0./sqrt(diag(Y))) )';
    V =   LC*Z*diag(1.0./sqrt(diag(Y)));

end

function y = assess(F,G,X,U,P,S,h,T,UN,VN,RN,UP,VP,RP)

    up = UP(1:RP,:);
    vp = VP(:,1:RP);

    p = up*P;

    un = UN(1:RN,:);
    vn = VN(:,1:RN);

    x = un*X;

    f = @(x,u,p) un*F(vn*x,u,vp*p);
    g = @(x,u,p)    G(vn*x,u,vp*p);

    y = ode(f,g,x,U,p,S,h,T);

end

function y = ode(f,g,x,u,p,S,h,E)

 T = (E-S)/h;
 y = zeros(numel(g(x,u(:,1),p)),T);

 for t=1:T
  x = x + h*f(x + (0.5*h)*f(x,u(:,t),p),u(:,t),p); %Improved Eulers Method
  y(:,t) = g(x,u(:,t),p);
 end

end

