function unified(s)
% unified (Version 1.2)
% by Christian Himpe, 2013 ( http://gramian.de )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%*

if(exist('emgr')~=2) disp('emgr framework is required. Download at http://gramian.de/emgr.m'); return; end

%%%%%%%% Setup %%%%%%%%

rand('seed',65537);

J = 10;
N = 100;
O = J;
R = O;
T = [0 0.01 1];
L = (T(3)-T(1))/T(2);

A = rand(N,N); A(1:N+1:end) = -0.55*N; A = 0.5*(A+A');
B = rand(N,J);
C = B';

U = [eye(J) zeros(J,L-J)];
X = J*ones(N,1);
P = rand(N,1)-0.5;
F = @(x,u,p) A*asinh(x)+B*u+p;
G = @(x,u,p) C*x;

 norm2 = @(y) sqrt(sum(sqrt(sum((y.*y),1)).^2));

y0 = rk2(F,G,T,X,U,P);

%%%%%%%% Balanced Truncation via Empirical Controllability Gramian and Empirical Observability Gramian %%%%%%%%

tic;
WC = emgr(F,G,[J N O],T,'c',P);
WO = emgr(F,G,[J N O],T,'o',P);
[UU D VV] = balance(WC,WO); UU = UU(1:R,:); VV = VV(:,1:R);
x = UU*X;
f = @(x,u,p) UU*F(VV*x,u,p);
g = @(x,u,p) G(VV*x,u,p);

y1 = rk2(f,g,T,x,U,P);
time_and_error_WCWO = [toc,norm2(y1-y0)./norm2(y0)]

%%%%%%%% Direct Truncation via Empirical Cross Gramian %%%%%%%%

tic;
WX = emgr(F,G,[J N O],T,'x',P);
[UU D VV] = svd(WX); UU = UU(:,1:R); VV = VV(:,1:R)';
x = VV*X;
f = @(x,u,p) VV*F(UU*x,u,p);
g = @(x,u,p) G(UU*x,u,p);

y2 = rk2(f,g,T,x,U,P);
time_and_error_WX = [toc,norm2(y2-y0)./norm2(y0)]

%%%%%%%% Parameter Reduction via Empirical Sensitivity Gramian %%%%%%%%

tic;
WS = emgr(F,G,[J N O],T,'s',P);
[PP D QQ] = svd(WS{2}); PP = PP(:,1:R); QQ = QQ(1:R,:);
p = QQ*P;
f = @(x,u,p) F(x,u,PP*p);
g = @(x,u,p) G(x,u,PP*p);
y3 = rk2(f,g,T,X,U,p);
time_and_error_WS = [toc,norm2(y3-y0)./norm2(y0)]

%%%%%%%% Parameter Reduction via Empirical Identifiability Gramian %%%%%%%%

tic;
WI = emgr(F,G,[J N O],T,'i',P);
[PP D QQ] = svd(WI{2}); PP = PP(:,1:R); QQ = QQ(:,1:R)';
p = QQ*P;
f = @(x,u,p) F(x,u,PP*p);
g = @(x,u,p) G(x,u,PP*p);
y4 = rk2(f,g,T,X,U,p);
time_and_error_WI = [toc,norm2(y4-y0)./norm2(y0)]

%%%%%%%% Combined Reduction (State and Parameter) via Empirical Joint Gramian %%%%%%%%

tic;
WJ = emgr(F,G,[J N O],T,'j',P);
[UU D VV] = svd(WJ{1}); UU = UU(:,1:R); VV = VV(:,1:R)';
[PP D QQ] = svd(WJ{2}); PP = PP(:,1:R); QQ = QQ(:,1:R)';
x = VV*X;
p = QQ*P;
f = @(x,u,p) VV*F(UU*x,u,PP*p);
g = @(x,u,p) G(UU*x,u,PP*p);
y5 = rk2(f,g,T,x,U,p);
time_and_error_WJ = [toc,norm2(y5-y0)./norm2(y0)]

%%%%%%%% Plot %%%%%%%%

if(nargin<1 || s==0 ) return; end

l = (1:-0.01:0)';
cmap = [l,l,ones(101,1)];	% blue-ish colormap

f1 = figure('units','normalized','position',[0 0 .35 .30]);
imagesc(abs(y0-y1)./abs(y0)); xlabel('Time'); ylabel('Relative Error'); colorbar; colormap(cmap); set(gca,'YDir','normal');
if(s==2) saveas(f1,'wcwo','epsc'); end

f2 = figure('units','normalized','position',[0 0 .35 .30]);
imagesc(abs(y0-y2)./abs(y0)); xlabel('Time'); ylabel('Relative Error'); colorbar; colormap(cmap); set(gca,'YDir','normal');
if(s==2) saveas(f2,'wx','epsc'); end

f3 = figure('units','normalized','position',[0 0 .35 .30]);
imagesc(abs(y0-y3)./abs(y0)); xlabel('Time'); ylabel('Relative Error'); colorbar; colormap(cmap); set(gca,'YDir','normal');
if(s==2) saveas(f3,'ws','epsc'); end

f4 = figure('units','normalized','position',[0 0 .35 .30]);
imagesc(abs(y0-y4)./abs(y0)); xlabel('Time'); ylabel('Relative Error'); colorbar; colormap(cmap); set(gca,'YDir','normal');
if(s==2) saveas(f4,'wi','epsc'); end

f5 = figure('units','normalized','position',[0 0 .35 .30]);
imagesc(abs(y0-y5)./abs(y0)); xlabel('Time'); ylabel('Relative Error'); colorbar; colormap(cmap); set(gca,'YDir','normal');
if(s==2) saveas(f5,'wj','epsc'); end

%%%%%%%% Integrator %%%%%%%%

function y = rk2(f,g,t,x,u,p)

 h = t(2);
 T = (t(3)-t(1))/h;
 y = zeros(numel(g(x,u(:,1),p)),T);

 for t=1:T
  x = x + h*f(x + 0.5*h*f(x,u(:,t),p),u(:,t),p); %Improved Eulers Method
  y(:,t) = g(x,u(:,t),p);
 end

%%%%%%%% Balancer %%%%%%%%

function [X Y Z] = balance(WC,WO)

 L = chol(WC+eye(size(WC,1)))-eye(size(WC,1));
 [U Y V] = svd(L*WO*L');
 X = diag(sqrt(diag(Y))) * V' / L';
 Z = L'*U*diag(1./sqrt(diag(Y)));
