function Z = acom(n)
% acom (Version 1.2)
% by Christian Himpe, 2013-2015 ( http://wwwmath.uni-muenster.de/u/himpe )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%*

format long;
warning off all;

if(nargin==0), n = 3; end;

d = 0.6;
M = 9;
Z = zeros(10,4,M);

%%%%%%%% GENERIC MODEL %%%%%%%%

rand('seed',1009);
randn('seed',1009);

J = n;
N = n*n;
O = n;

A = rand(N,N);
B = rand(N,J);
C = rand(O,N);
D = 0;
F = 0;
P = A(:); P(1:N+1:end) = [];

EN = speye(N*N); EX = ~eye(N); EX = EX(:); EN(1:N*N+1:end) = EN(1:N*N+1:end).*EX'; EN(:,all(~EN,1) ) = [];
ED = -0.5*N*speye(N);

T = [0,0.2,20.0];
L = (T(3)-T(1))/T(2);
X = zeros(N,1);
U = [ones(J,1),zeros(J,L-1)];

Q = (1/N)*ones(N*N-N,1);
S = 1.0;
K = @(p) reshape(EN*p,[N,N])+ED;

Y = ode(K(P),B,C,D,F,T,X,U);
Y = Y + 0.01*max(max(abs(Y)))*randn(O,L);
%figure; plot(1:L,Y);

%
Z(1,:,:) = redopt(K,B,C,D,F,X,U,Q,T,S,Y,P,-1,M,0.01); fprintf('#');                     % FULL
Z(2,:,:) = redopt(K,B,C,D,F,X,U,Q,T,S,Y,P,[0,2,0,0,0,0.02,0,0,0,1,0,0],M,0.01); fprintf('#');  % ORIGINAL
Z(3,:,:) = redopt(K,B,C,D,F,X,U,Q,T,S,Y,P,[0,2,0,0,0,0.02,0.02,1,0,1,0,0],M,0.01); fprintf('#');  % DATA-DRIVEN
Z(4,:,:) = redopt(K,B,C,D,F,X,U,Q,T,S,Y,P,[0,2,0,0,0,0.02,0,0,1,1,0,0],M,0.01); fprintf('#');  % MONTE-CARLO
Z(5,:,:) = redopt(K,B,C,D,F,X,U,Q,T,S,Y,P,[0,2,0,0,0,0.02,0.01,1,1,1,0,0],M,0.01); fprintf('#'); % MONTE-CARLO+DATA-DRIVEN
%

fprintf('|');

%%%%%%%% CONNECTIVITY MODEL %%%%%%%%

rand('seed',1009);
randn('seed',1009);

J = n;
N = n*n*5;
O = n*n;

B = sparse(N,J); B(1:n*n,1:J) = rand(n*n,J);
C = sparse(O,N);

for m=0:O-1
	C(m+1,n*n+4*m+3) = 1.0;
	C(m+1,n*n+4*m+4) = -1.0;
end

X = zeros(N,1);
U = [ones(J,1),zeros(J,L-1)];

A = effconn_1(O);
K = @(p) effconn_2(A,O,p,EN,ED);

Y = ode(K(P),B,C,D,F,T,X,U);
Y = Y + 0.01*max(max(abs(Y)))*randn(O,L);
%figure; plot(1:L,Y);

%
Z(6,:,:)  = redopt(K,B,C,D,F,X,U,Q,T,S,Y,P,-1,M,0.0001); fprintf('#');                    % FULL
Z(7,:,:)  = redopt(K,B,C,D,F,X,U,Q,T,S,Y,P,[0,2,2,0,5,0.000068,0.0,0,0,1,0,0],M,0.0001); fprintf('#'); % ORIGINAL
Z(8,:,:)  = redopt(K,B,C,D,F,X,U,Q,T,S,Y,P,[0,2,2,0,5,0.000068,0.000074,1,0,1,0,0],M,0.0001); fprintf('#'); % DATA-DRIVEN
Z(9,:,:)  = redopt(K,B,C,D,F,X,U,Q,T,S,Y,P,[0,2,2,0,5,0.000068,0.0,0,1,1,0,0],M,0.0001); fprintf('#'); % MONTE-CARLO
Z(10,:,:) = redopt(K,B,C,D,F,X,U,Q,T,S,Y,P,[0,2,2,0,5,0.000068,0.000068,1,1,1,0,0],M,0.0001); fprintf('#'); % MONTE-CARLO+DATA-DRIVEN
%

save('acom.mat','Z','M');

end

%%%%%%%% OPTIMIZER %%%%%%%%
function z = redopt(K,B,C,D,F,X,U,Q,T,S,Y,P,nf,n,reg)

	norm2   = @(y) sqrt(sum(sum(y.*y)));
	norm2T  = @(y) sqrt(T(2)*sum(sum(y.*y)));
	norm22T = @(y) T(2)*sum(sum(y.*y));

	z = zeros(1,4,n);

	if(nf==-1),
		z(1,1,:) = NaN;

		tic;
		j = @(p) norm22T(Y-ode(K(p),B,C,D,F,T,X,U))+reg*(p'*p);
		p = fminunc(j,Q,optimset('Display','off'));
		z(1,2,:) = toc;

		y = ode(K(p),B,C,D,F,T,X,U);

		z(1,3,:) = norm2T(Y - y)/norm2T(Y);
		z(1,4,:) = norm2(P - p)/norm2(P);
	else,
		global gt;
		global gd;
		gt = tic;

		gd = zeros(1,1,n);
		XP = optmor(K,B,C,D,F,T,n,X,U,Q,S,nf,Y);
		z(1,1,:) = gd;

		for I=2:n,
			XX = XP{I,1};
			PP = XP{I,2};
			q = PP'*Q;
			b = XX'*B;
			f = XX'*F;
			c = C*XX;
			x = XX'*X;

			tic;
			j = @(p) norm22T(Y-ode(XX'*K(PP*p)*XX,b,c,0,0,T,x,U))+reg*( (p'*PP')*(PP*p) );
			p = fminunc(j,q,optimset('Display','off'));
			z(1,2,I) = toc;

			y = ode(XX'*K(PP*p)*XX,b,c,0,0,T,x,U);

			z(1,3,I) = norm2T(Y - y)/norm2T(Y);
			z(1,4,I) = norm2(P - PP*p)/norm2(P);
		end;
	end;
end

%%%%%%%% TRANSFORM %%%%%%%%
function a = effconn_1(N)

	tau_s = 0.65;
	tau_f = 0.41;
	tau_0 = 0.98;
	  E_0 = 0.34;
	alpha = 0.32;

	H = [-tau_s,-tau_f,0,0;...
             1.0,0,0,0;...
	     0,(1.0-(1.0-E_0)*(1.0-log(1.0-E_0)))/(tau_0*E_0),-1.0/tau_0,(1.0-alpha)/(tau_0*alpha);...
	     0,(1.0/tau_0),0,-1.0/(tau_0*alpha)];

	a = sparse(5*N,5*N);

	for M=0:N-1,
		t = N+4*M+1;
		a(t:t+3,t:t+3) = H;
		a(t,M+1) = 1.0;
	end;
end

function A = effconn_2(A,N,p,EN,ED)

	A(1:N,1:N) = reshape(EN*p,[N,N])+ED;
end

%%%%%%%% INTEGRATOR %%%%%%%%
function y = ode(A,B,C,D,F,T,X,U)

	h = T(2);
	H = 1.0./h;
	L = round((T(3)-T(1))/h);

	if(exist('OCTAVE_VERSION'))
		f = @(y,t) A*y + B*U(:,1.0+min(round(t*H),L-1)) + F;
		y = C*lsode(f,X,linspace(0,h*L,L))'; % + D*U;
	else
		f = @(t,y) A*y + B*U(:,1.0+min(round(t*H),L-1)) + F;
		y = C*deval( ode45(f,[0,h*L],X), linspace(0,h*L,L) ); % + D*U;
	end
end
