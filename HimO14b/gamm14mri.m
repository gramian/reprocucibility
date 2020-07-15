function Z = gamm14mri(n)
% gamm14mri (Version 1.0)
% by Christian Himpe, 2014 ( http://wwwmath.uni-muenster.de/u/himpe )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%*

rand('seed',1009);

%%%%%%%% MRI %%%%%%%%

%%%% Setup %%%%

N = 5*n;
M = n;
O = n;
P = n*n;

p = rand(n);
p(1:n+1:end) = -0.5*n;
p = p(:);

S = 0;
h = 0.01;
T = 20.0;

 A0 = mri_dyn_lin_0(n);
 A = @(p) A0 + mri_dyn_lin_p(p);
 B = [eye(n);zeros(N-n,n)];
 C = [zeros(n,n+n+n),eye(n),-eye(n)];

X = [zeros(2*n,1);ones(3*n,1)];
u = @(t) ones(n,1)*(t<=h);

%%%% Original %%%%

 tic;
 Y = integrate(@(x,u,p) A(p)*x+B*u,@(x,u,p) C*x,S,h,T,X,u,p);
 ON_MRI_ORIG = toc

%%%% Reduction %%%%

%% Gramian-Based %%

 tic;
 WJ = emgr(@(x,u,p) A(p)*x+B*u,@(x,u,p) C*x,[M N O],[S h T],'j',p,[0 0 3 3 3 3 0 0 0 -1 0 0],1,0,X);
 [GU D GV] = svd(WJ{1});
 [GP D GQ] = svd(WJ{2});
 OFF_MRI_EMGR = toc

%% Optimization-Based %%

 tic;
 PV = optmor(p,A,B,C,0,0,[S h T],max(N,P),X,1,1,[0,3,1,2,2,0,-1,0,0],Y);
 OP = PV{1};
 OV = PV{2};
 OFF_MRI_OPTMOR = toc

%%%% Assessment %%%%

 mri_l2 = zeros(P,N,2);
 mri_h2 = zeros(P,N,2);

 for I=1:P
 for J=1:N
     TU = GU(:,1:J);
     TV = TU';
     TP = GP(:,1:I); 
     TQ = TP';

     x = TV*X;
     q = TQ*p;
     a = TV*A(TP*q)*TU;
     b = TV*B;
     c = C*TU;
     mri_dyn_red = @(x,u,p) a*x+b*u;
     mri_out_red = @(x,u,p) c*x;

     y = integrate(mri_dyn_red,mri_out_red,S,h,T,x,u,q);

     mri_l2(I,J,1) = sqrt(sum(sum((Y-y).*(Y-y))))./sqrt(sum(sum(Y.*Y)));

     TU = OV(:,1:J);
     TV = TU';
     TP = OP(:,1:I);
     TQ = TP';

     x = TV*X;
     q = TQ*p;
     a = TV*A(TP*q)*TU;
     b = TV*B;
     c = C*TU;
     mri_dyn_red = @(x,u,p) a*x+b*u;
     mri_out_red = @(x,u,p) c*x;

     y = integrate(mri_dyn_red,mri_out_red,S,h,T,x,u,q);

     mri_l2(I,J,2) = sqrt(sum(sum((Y-y).*(Y-y))))./sqrt(sum(sum(Y.*Y)));
 end
 end

%%%% Plot %%%%

 mri_l2 = min(mri_l2,1.0);

 figure;

 cmap = colormap('autumn'); cmap = flipud(cmap);

 subplot(1,2,1); 
 h = surf(mri_l2(:,:,1)); shading interp; view(150,30);
 colormap(cmap); set(gca,'clim',[0 0.2]);  set(h','edgecolor','k');
 xlim([1 N]); ylim([1 P]); zlim([0.001 1]); set(gca,'zscale','log');
 title('Relative L2 Output Error (emgr)','fontweight','bold','fontsize',12);
 xlabel('States','fontsize',12); 
 ylabel('Parameters','fontsize',12);

 subplot(1,2,2);
 h = surf(mri_l2(:,:,2)); shading interp; view(150,30);
 colormap(cmap); set(gca,'clim',[0 0.2]); set(h','edgecolor','k');
 xlim([1 N]); ylim([1 P]); zlim([0.001 1]); set(gca,'zscale','log');
 title('Relative L2 Output Error (optmor)','fontweight','bold','fontsize',12); 
 xlabel('States','fontsize',12); 
 ylabel('Parameters','fontsize',12);

print -depsc 'mri2.eps'
end


%%%%%%%% Local %%%%%%%%

function A = mri_dyn_lin_p(p)

 n = sqrt(numel(p));

 A = sparse(5*n,5*n);
 A(1:n,1:n) = reshape(p,[n,n]);

end

function A = mri_dyn_lin_0(n)

 k = 0.65;
 g = 0.41;
 a = 0.32;
 r = 0.34;
 t = 0.98;

 A = sparse(5*n,5*n);

 A(n+1:n+n,1:n) = speye(n); 
 A(n+1:n+n,n+1:n+n) = -k*speye(n);
 A(n+1:n+n,n+n+1:n+n+n) = -g*speye(n);

 A(n+n+1:n+n+n,n+1:n+n) = speye(n);

 A(n+n+n+1:n+n+n+n,n+n+1:n+n+n) = (1.0/t)*speye(n);
 A(n+n+n+1:n+n+n+n,n+n+n+1:n+n+n+n) = (-1.0/(t*a))*speye(n);

 A(n+n+n+n+1:n+n+n+n+n,n+n+1:n+n+n) = (1.0/(t*r))*(1.0-(1.0-r)*(1.0-log(1.0-r)))*speye(n);
 A(n+n+n+n+1:n+n+n+n+n,n+n+n+1:n+n+n+n) = (-(1.0-a)/(a*t))*speye(n);
 A(n+n+n+n+1:n+n+n+n+n,n+n+n+n+1:n+n+n+n+n) = (-1.0/t)*speye(n);

end

function y = integrate(f,g,S,h,T,x,u,p)

if(exist('OCTAVE_VERSION')),
    x = lsode(@(y,t) f(y,u(t),p),x,linspace(S,T,T/h))';
else,
    x = deval(ode45(@(t,y) f(y,u(t),p),[S T],x),linspace(S,T,T/h));
end;

y = g(x,0,p);

end

