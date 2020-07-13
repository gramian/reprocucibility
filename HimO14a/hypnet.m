function z = hypnet()
% hyperbolic network reduction
% by Christian Himpe, 2013-2014 ( http://gramian.de )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%*

rand('seed',1009);
if(exist('emgr')~=2) disp('emgr framework is required. Download at http://gramian.de/emgr.m'); return; end

%%%%%%%% Setup %%%%%%%%

 J = 8;
 O = J;
 N = 64;
 T = [0 0.01 1.0];
 L = (T(3)-T(1))/T(2);
 U = [ones(J,1) zeros(J,L-1)];
 X =  zeros(N,1);

 %A = network(N,2);
 B = rand(N,J);
 C = rand(O,N);

 P = zeros( ((N*N)-N)/2, 1 );

 LIN = @(x,u,p) tvs(L,N)*x + B*u;
 ELN = @(x,u,p) sps(p,N)*x + [C' B]*u;
 OUT = @(x,u,p) C*x;
 EOU = @(x,u,p) [C;B']*x;

 norm2 = @(y) sqrt(sum(sum(y.*y)));

%%%%%%%% Reduction %%%%%%%%

% FULL
 rand('seed',1009);
 tic; Y = rk2(LIN,OUT,T,X,U,P); z(1) = toc;

% OFFLINE
 tic;
 WJ = emgr(ELN,EOU,[J+O N O+J],T,'j',P,0,1,0,1);
 [UU D VV] = svd(WJ{1}); UU = UU(:,1:19); VV = VV(:,1:19)'; %diag(D)'
 [PP E QQ] = svd(WJ{2}); PP = PP(:,1:65); QQ = QQ(:,1:65)'; %diag(E)'
 b = VV*B;
 c = C*UU;
 x = VV*X;
 p = QQ*P;
 lin = @(x,u,p) VV*(tvs(L,N)*(UU*x)) + b*u;
 out = @(x,u,p) c*x;
 z(2) = toc;

% ONLINE
 rand('seed',1009);
 tic; y = rk2(lin,out,T,x,U,PP*p); z(3) = toc;

%%%%%%%% Output %%%%%%%%

% TERMINAL
 z(4) = norm2(norm2(Y - y)./norm2(Y));
 RELER = abs(Y - y)./abs(Y);

% PLOT
 l = (1:-0.01:0)'; cmap = [l,l,ones(101,1)];
 figure;
 semilogy(diag(D),'LineWidth',2); xlim([1 N]); xlabel('Singular Value'); ylabel('State');
 print -depsc cross.epsc;

 figure;
 semilogy(diag(E),'LineWidth',2); xlim([1 (N*N-N)/2]); xlabel('Singular Value'); ylabel('Parameter');
 print -depsc ident.epsc;

 figure;
 imagesc(Y); caxis([0 max(max(Y))]); colorbar; colormap(cmap); set(gca,'YTick',1:N); xlabel('Time'); ylabel('States');
 print -depsc hypfull.epsc;

 figure;
 imagesc(y); caxis([0 max(max(Y))]); colorbar; colormap(cmap); set(gca,'YTick',1:N); xlabel('Time'); ylabel('States');
 print -depsc hypred.epsc;

 figure;
 imagesc(RELER); caxis([0 max(max(RELER))]); colorbar; colormap(cmap); set(gca,'YTick',1:N); xlabel('Time'); ylabel('Relative Error');
 print -depsc hypnet.epsc;

%%%%%%%% Parametrized Symmetric System %%%%%%%%

function A = sps(p,N)

 persistent LL; if(isempty(LL)), LL = tril(ones(N),-1); end;

 A = -N*eye(N);
 A(LL==1) = p;
 A = 0.01*(A+A');  

%%%%%%%% Time Varying Symmetric System %%%%%%%%

function A = tvs(L,N)

 persistent TI; if(isempty(TI)), TI = 1; end;
 persistent PR; if(isempty(PR)), PR = zeros(N,N); end;
 persistent NL; if(isempty(NL)), NL = zeros(N,2); end;
 
 if(TI<=N)
  NL(TI,:) = [2.0*log(TI/2.0),2.0*pi*rand(1)];

  for I=1:(TI-1)
   S = pi-abs(pi-abs(NL(TI,2)-NL(I,2)));
   PR(TI,I) = ( NL(I,1) + 2.0*log(S*0.5) ) < 0;
  end
 end

 A = -N*eye(N);
 A = A + PR;
 A = 0.01*(A+A');

 TI = TI + 1;
 if(TI==L+1), TI = 1; PR = zeros(N,N); NL = zeros(N,2); end

%%%%%%%% Network %%%%%%%%

function A = network(N,v)

 d = 1;
 b = 1;

 a = zeros(N,1+d);
 A = zeros(N,N);

 for n=2:N
  r = 2*b*log(n/v);
  t = 2*pi*rand(1,d);
  a(n,:) = [r,t];

  for m=1:n-1
   s = pi-norm(pi-norm(t-a(m,2:d+1)));
   A(n,m) = ( a(m,1) + 2*b*log(s/2) ) < 0;
  end
 end

 A = 0.5*(A+A');

%%%%%%%% Integrator %%%%%%%%

function y = rk2(f,g,t,x,u,p)

 h = t(2);
 T = (t(3)-t(1))/h;
 y = zeros(numel(g(x,u(:,1),p)),T);

 for t=1:T
  x = x + h*f(x + 0.5*h*f(x,u(:,t),p),u(:,t),p); %Improved Eulers Method
  y(:,t) = g(x,u(:,t),p);
 end
 
