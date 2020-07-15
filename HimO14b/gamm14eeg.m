function Z = gamm14eeg(n)
% gamm14eeg (Version 1.0)
% by Christian Himpe, 2014 ( http://wwwmath.uni-muenster.de/u/himpe )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%*

rand('seed',1009);

%%%%%%%% EEG %%%%%%%%

%%%% Setup %%%%

N = 10*n;
M = n;
O = n;
P = 3*n*n;

af = rand(n); af(1:n+1:end) = 0;
ab = rand(n); ab(1:n+1:end) = 0;
al = rand(n); al(1:n+1:end) = 0;
p = [af(:);ab(:);al(:)];

S = 0;
h = 0.001;
T = 2.0;

 A0 = eeg_dyn_lin_0(n);
 A = @(p) A0 + eeg_dyn_lin_p(p);
 B = [zeros(7*n,n);8*eye(n);zeros(2*n,n)];
 C = [eye(n),-eye(n),zeros(n,8*n)];

X = zeros(N,1);
u = @(t) ones(n,1)*(t<=h);

%%%% Original %%%%

 tic;
 Y = integrate(@(x,u,p) A(p)*x+B*u,@(x,u,p) C*x,S,h,T,X,u,p);
 ON_EEG_ORIGINAL = toc

%%%% Reduction %%%%

%% Gramian-Based %%

 tic;
 WJ = emgr(@(x,u,p) A(p)*x+B*u,@(x,u,p) C*x,[M N O],[S h T],'j',p,[1 0 3 3 3 3 0 0 0 -1 0 0],u,0,X);
 [GU D GV] = svd(WJ{1});
 [GP D GQ] = svd(WJ{2});
 OFF_EEG_EMGR = toc

%% Optimization-Based %%

 tic;
 PV = optmor(p,A,B,C,0,0,[S h T],max(N,P),X,1,1,[0,1,1,2,2,0,-1,0,0],Y);
 OP = PV{1};
 OV = PV{2};
 OFF_EEG_OPTMOR = toc

%%%% Assessment %%%%

 eeg_l2 = zeros(P,N,2);
 eeg_h2 = zeros(P,N,2);

 for I=1:P
 for J=1:N
     TU = GU(:,1:J);
     TV = TU';
     TP = GP(:,1:I); 
     TQ = TP';

     x = TV*zeros(N,1);
     q = TQ*p;
     a = TV*A(TP*q)*TU;
     b = TV*B;
     c = C*TU;
     eeg_dyn_red = @(x,u,p) a*x+b*u;
     eeg_out_red = @(x,u,p) c*x;

     y = integrate(eeg_dyn_red,eeg_out_red,S,h,T,x,u,q);

     eeg_l2(I,J,1) = sqrt(sum(sum((Y-y).*(Y-y))))./sqrt(sum(sum(Y.*Y)));

     TU = OV(:,1:J);
     TV = TU';
     TP = OP(:,1:I);
     TQ = TP';

     x = TV*zeros(N,1);
     q = TQ*p;
     a = TV*A(TP*q)*TU;
     b = TV*B;
     c = C*TU;
     eeg_dyn_red = @(x,u,p) a*x+b*u;
     eeg_out_red = @(x,u,p) c*x;

     y = integrate(eeg_dyn_red,eeg_out_red,S,h,T,x,u,q);

     eeg_l2(I,J,2) = sqrt(sum(sum((Y-y).*(Y-y))))./sqrt(sum(sum(Y.*Y)));
 end
 end

%%%% Plot %%%%

 eeg_l2 = min(eeg_l2,1.0);

 figure;

 cmap = colormap('autumn'); cmap = flipud(cmap);

 subplot(1,2,1); 
 h = surf(eeg_l2(:,:,1)); shading interp; view(150,30);
 colormap(cmap); set(gca,'clim',[0 0.2]);  set(h','edgecolor','k');
 xlim([1 N]); ylim([1 P]); zlim([0.001 1]); set(gca,'zscale','log');
 title('Relative L2 Output Error (emgr)','fontweight','bold','fontsize',12);
 xlabel('States','fontsize',12); 
 ylabel('Parameters','fontsize',12);

 subplot(1,2,2);
 h = surf(eeg_l2(:,:,2)); shading interp; view(150,30);
 colormap(cmap); set(gca,'clim',[0 0.2]); set(h','edgecolor','k');
 xlim([1 N]); ylim([1 P]); zlim([0.001 1]); set(gca,'zscale','log');
 title('Relative L2 Output Error (optmor)','fontweight','bold','fontsize',12); 
 xlabel('States','fontsize',12); 
 ylabel('Parameters','fontsize',12);

 print -depsc 'eeg2.eps'
end

%%%%%%%% Local %%%%%%%%

function A = eeg_dyn_lin_p(p)

 m = sqrt(numel(p)/3);
 n = 5*m;
 N = 10*m;
 M = m*m;

KE = 1/8.0;  HE = 4.0;  KEHE = KE*HE;
KI = 1/16.0; HI = 32.0; KIHI = KI*HI;

SIG = 0.14;

 AF = reshape(p(1:M),[m,m]);
 AB = reshape(p(M+1:M+M),[m,m]);
 AL = reshape(p(M+M+1:end),[m,m]);

 A = sparse(N,N);

A(n+1:n+m,2*m+1:3*m) = KEHE*SIG*(AB+AL);

A(n+m+m+1:n+m+m+m,1:m) = KEHE*SIG*(AF+AL);
A(n+m+m+1:n+m+m+m,m+1:m+m) = -KEHE*SIG*(AF+AL);

A(n+m+m+m+1:n+m+m+m+m,1:m) = KEHE*SIG*(AB+AL);
A(n+m+m+m+1:n+m+m+m+m,m+1:m+m) = -KEHE*SIG*(AB+AL);

end

function A = eeg_dyn_lin_0(m)

n = 5*m;
N = 10*m;

KE = 1/8.0;  HE = 4.0;  KE2 = 2.0*KE; KEKE = KE*KE; KEHE = KE*HE;
KI = 1/16.0; HI = 32.0; KI2 = 2.0*KI; KIKI = KI*KI; KIHI = KI*HI;

G1 = 1.0;
G2 = 1.0;
G3 = 0.50;
G4 = 0.50;
G5 = 0.25;

SIG = 0.14;

A = sparse(N,N);

A(1:n,n+1:N) = speye(n);

A(n+1:n+m,2*m+1:3*m) = KEHE*G2*SIG*speye(m);
A(n+1:n+m,5*m+1:6*m) = -KE2*speye(m);
A(n+1:n+m,0*m+1:1*m) = -KEKE*speye(m);

A(n+m+1:n+m+m,3*m+1:4*m) = KIHI*G4*SIG*speye(m);
A(n+m+1:n+m+m,4*m+1:5*m) = -KIHI*G4*SIG*speye(m);
A(n+m+1:n+m+m,6*m+1:7*m) = -KI2*speye(m);
A(n+m+1:n+m+m,1*m+1:2*m) = -KIKI*speye(m);

A(n+m+m+1:n+m+m+m,0*m+1:1*m) = KEHE*G1*SIG*speye(m);
A(n+m+m+1:n+m+m+m,1*m+1:2*m) = -KEHE*G1*SIG*speye(m);
A(n+m+m+1:n+m+m+m,7*m+1:8*m) = -KE2*speye(m);
A(n+m+m+1:n+m+m+m,2*m+1:3*m) = -KEKE*speye(m);

A(n+m+m+m+1:n+m+m+m+m,0*m+1:1*m) = KEHE*G3*SIG*speye(m);
A(n+m+m+m+1:n+m+m+m+m,1*m+1:2*m) = -KEHE*G3*SIG*speye(m);
A(n+m+m+m+1:n+m+m+m+m,8*m+1:9*m) = -KE2*speye(m);
A(n+m+m+m+1:n+m+m+m+m,3*m+1:4*m) = -KEKE*speye(m);

A(n+m+m+m+m+1:n+m+m+m+m+m,3*m+1:4*m) = KIHI*G5*SIG*speye(m);
A(n+m+m+m+m+1:n+m+m+m+m+m,4*m+1:5*m) = -KIHI*G5*SIG*speye(m);
A(n+m+m+m+m+1:n+m+m+m+m+m,9*m+1:10*m) = -KI2*speye(m);
A(n+m+m+m+m+1:n+m+m+m+m+m,4*m+1:5*m) = -KIKI*speye(m);

end

function y = integrate(f,g,S,h,T,x,u,p)

if(exist('OCTAVE_VERSION')),
    x = lsode(@(y,t) f(y,u(t),p),x,linspace(S,T,T/h))';
else,
    x = deval(ode45(@(t,y) f(y,u(t),p),[S T],x),linspace(S,T,T/h));
end;

y = g(x,0,p);

end

