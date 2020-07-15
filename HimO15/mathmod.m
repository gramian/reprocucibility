function mathmod(N)
% mathmod 2015 (Companion Code)
% by Christian Himpe, 2014 ( http://gramian.de )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%*

if(exist('emgr')~=2) disp('emgr framework is required. Download at http://gramian.de/emgr.m'); return; end

rand('seed',1009);
randn('seed',1009);

%%%%%%%% Setup %%%%%%%%

J = 1;
if(nargin<1), N = 1000; end;
O = 1;

S = 0;
h = 0.01;
T = 1.0;

gn = @(x,y) exp(10*x) + y - 1.0;
A1 = spdiags(ones(N-1,1),-1,N,N) - speye(N);
A2 = spdiags([ones(N-1,1);0],0,N,N) - spdiags(ones(N,1),1,N,N);
Hn = @(x,u,p) gn(A1*x,A1*(p.*x)) - gn(A2*x,A2*(p.*x)) + [u;sparse(N-1,1)];
Gn = @(x,u,p) x(1);
Xn = zeros(N,1);
Un = exp(-10.0*[h:h:T]);
PL = repmat(linspace(0.5,1.5,10),[N,1]);
Pl = rand(N,1)+0.5;

l2norm = @(m) sqrt(sum(sum(m.*m,1)));
l8norm = @(m) max(sqrt(sum(m.*m,1)));

ONLINE  = zeros(1,3);
OFFLINE = zeros(1,2);
SINGVAL = zeros(2,N);
L2ERROR = zeros(2,N-1);
L8ERROR = zeros(2,N-1);

tic;
Y = assess(Hn,Gn,Xn,Un,Pl,S,h,T,speye(N),speye(N),N,speye(N),speye(N),N);
ONLINE(3) = toc;

%%%%%%%% Nonlinear Benchmark State Reduction %%%%%%%%

% BT
tic;
WC = emgr(Hn,Gn,[1 N 1],[S h T],'c',PL,[0,0,0,0,0,0,0,0,0,0,0,1],Un,[],[],[],0.01);
WO = emgr(Hn,Gn,[1 N 1],[S h T],'o',PL,[0,0,0,0,0,0,0,0,0,0,0,1],Un,[],[],[],0.01);
[UB DB VB] = bal(WC,WO);
OFFLINE(1) = toc;
%SINGVLA(1,:) = diag(DB);
ONLINE(1) = 0;
for I=1:N-1, tic; y = assess(Hn,Gn,Xn,Un,Pl,S,h,T,UB,VB,I,speye(N),speye(N),N); ONLINE(1) = ONLINE(1) + toc; L2ERROR(1,I) = l2norm(Y-y); L8ERROR(1,I) = l8norm(Y-y); end;
ONLINE(1) = ONLINE(1)/(N-1);
fprintf('#');

% WX
tic;
WX = emgr(Hn,Gn,[1 N 1],[S h T],'x',PL,[0,0,0,0,0,0,0,0,0,0,0,1],Un,[],[],[],0.01);
[UX DX VX] = svd(WX);
OFFLINE(2) = toc;
%SINGVAL(2,:) = diag(DX);
ONLINE(2) = 0;
for I=1:N-1, tic; y = assess(Hn,Gn,Xn,Un,Pl,S,h,T,UX',UX,I,speye(N),speye(N),N); ONLINE(2) = ONLINE(2) + toc; L2ERROR(2,I) = l2norm(Y-y); L8ERROR(2,I) = l8norm(Y-y); end;
ONLINE(2) = ONLINE(2)/(N-1);
fprintf('#\n');

OFFLINE

%%%% %%%% %%%% %%%%

ORIGL2 = l2norm(Y);
ORIGL8 = l8norm(Y);

save('mathmod.mat','N','OFFLINE','ONLINE','SINGVAL','L2ERROR','L8ERROR','ORIGL2','ORIGL8');

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

