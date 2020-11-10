
% e = [-3i-3 -1i+1 2i-2];
% e = [-10 -1 -0.01-i 10-3i  i 5];
% e = [-1 -2 3 7 6]
% % 
% A = diag(e);
% path = '../../results/projection_stability_ellipsoid/4by4_diag_real/'
% A(2,1) = -5;
%  A(1,2) = -5;


% %A is hermitian, complex and no instable eigenvalue associated to off-diagonal
% %elements
% e = [-10 -1 -6 10 5];
% A = diag(e);
% % A(1,3) = 1-i;
% % A(1,2) = -i;
% % A(2,1) = i;
% % A(3,1) = 1+i;
% % A(3,2) = -1-2i;
% % A(2,3) = -1+2i;
% A(1,3) = 10;
% A(1,2) = 3;
% A(2,1) = 3;
% A(2,3) = -1;
% A(1,4) = -10;
% A(1,5) = 5;
% A(2,4) = -10;
% % A(2,5) = 1;
% A(3,4) = -1;
% % A(3,5) = 3;
% % A(4,5) = -7;

% % 3 by 3 symmetric matrix
% e = [-2 -1];
% A = diag(e);
% A(1,2)=-i;
% A= (A + A')/2

% % 10 by 10 matrix
% e = [-10 -3 -6 0.1 0.1 0.1 0.1 0.1 -2 -4];
% A = diag(e);
% A(1,2) = 3;
% A(1,3) = 10;
% A(1,4) = -10;
% A(1,5) = 5;
% A(1,6) = -10;
% A(1,7) = 5;
% A(1,8) = -10;
% A(1,9) = 5;
% A(1,10) = -10;
% A(2,3) = -1;
% A(2,4) = -10;
% A(2,5) = 1;
% A(2,6) = 3;
% A(2,7) = 5;
% A(2,8) = -6;
% A(2,9) = -7;
% A(2,10) =-8;
% A(3,4) = -10;
% A(3,5) = -9;
% A(3,6) = 13;
% A(3,7) = 0;
% A(3,8) = 0;
% A(3,9) = 5;
% A(3,10) =5;
% A(4,5) = 1;
% A(4,6) = 1;
% A(4,7) = 1;
% A(4,8) = 1;
% A(4,9) = 5;
% A(4,10) =5;
% A(5,6) = -2;
% A(5,7) = -3;
% A(5,8) = -4;
% A(5,9) = -5;
% A(5,10) =9;
% A(6,7) = -1;
% A(6,8) = 1;
% A(6,9) = -3;
% A(6,10) =8;
% 
% A = A + A';
% path = '../../results/projection_stability_ellipsoid/10by10_sym_realroots/sqrt(n)/'


% A is hermitian, complex and one instable eigenvalue associated to off-diagonal
% % elements (Stable-draft : CASE1)
% e = [i 1-i 1 0 -1];
% A = diag(e);
% A(1,3) = -i;
% A(1,2) = -i;
% A(1,4) = i;
% A(1,5) = -i;
% A(3,1) = 1;
% A(4,2) = 1;
% A(3,2) = -1;
% A = (A + A')/2
% % % path = '../../results/projection_stability_ellipsoid/5by5_hermitian/'

% % A is hermitian, complex and one instable eigenvalue associated to off-diagonal
% % elements (Stable-draft : CASE1)
% e = [i 1-i 1];
% A = diag(e);
% A(1,3) = -i;
% A(1,2) = -i;
% A(3,1) = 1;
% A(3,2) = -1;
% % path = 'D:\Documents and Settings\fxorban\My Documents\fxorban\notes\Articles\figures\'

% %error 1
% A = [1 3 ;1 1]
% path = '../../results/projection_stability_ellipsoid/2by2_nonnormal/'

%error 2
% A = [1 3 2;1 1 4; 0 3 8]
% A = [1 8 7; 9 4 2; 0 3 8];
% path = '../../results/projection_stability_ellipsoid/3by3_nonnormal_methodedirecte/'

% A = [2 -1;1 2] % normal matrix, it worked!

% %next example results in "3by3_normal_cmplxroots"
% L = diag([0.5+i -1+1i/2 -2/3-1i/3 1/3-1i/2 1-1i/3]);
% W = [4 1+1i -1 3i 2-5i;1-1i 1 -1i 0.5 0.23;-1 1i 5 0.4-2i 1; -3i 0.5 0.4+2i 1 0; 2+5i 0.23 1 0 0.5];%is positive semi-definite
% [U,Sv,V] = svd(W);
% A = U*L*U';
% path = '../../results/projection_stability_ellipsoid/Discrete/5by5_normal_cmplxroots/'

%next example results in "6by6_normal_cmplxroots"
%Rate of Convergence (roc) is not appliable to this examples as 
%the eigenvalues saved in stock.mat permutes their column at some
%iterations.
% e = ones(6,1);
% A = full(spdiags(e,1,6,6));
% A(6,1) = 1;
% A = A*10000000;
% path = '../../results/projection_stability_ellipsoid/6by6_normal_cmplxroots/'

% %next example enabled me to correct a mistake; indeed, at one point lambda
% %gets negative when the starting points are closer to the stability limit
% %than all the eigenvalues of the matrix; => I had to prevent lambda from
% %being negative in find_lambda.m
% A = [-10 1  1  1;...
%       1 -5  1  1;...
%       1  1 -10 1;...
%       1  1  1 -10];


% A is anti-hermitian
% A = zeros(3,3);
% A(1,3) = 1;
% A(1,2) = -1;
% A(2,3) = -3;
% A = A - A';
% A(1,1) = -10;
% A(2,2 ) = 10;

% %Stable Matrix (should be found in one iteration)
% A = [-2+3i 0 0;0 -4 0;0 0 -1-2i]

%Heighth case (continuous)
delta = -0.1;
n =4;
A = diag(ones(n-1,1),-1);
A(1,n) = delta;
% s1 = (0.1)^(1/4) * exp(1i*3*pi/4);
% s2 = (0.1)^(1/4) * exp(1i*-3*pi/4);
% s = poly([s1 s1 s2 s2]);
% S = diag(ones(n-1,1),-1);
% S(:,n) = -flipud(s(2:n+1).');
% path = 'D:\Documents and Settings\fxorban\My Documents\Dropbox\Doctorat\Thesis\figures'

% %Generating random stable matrices(continuous)
% s = 3785 %3785 %2456
% rand('twister',s)
% si = 5;
% A = (rand(si)-0.7)*10;
% maxl = abs(max(real(eig(A))));
% S = A - 1.2*maxl*eye(si)-i*rand(si)
% A = A - 0.5*maxl*eye(si)
% % path = '../../results/projection_stability_ellipsoid/Continuous/6by6_randomized_real_IPAM2010/'

% %Generating random stable matrices (discrete)
% s = 3785 %3785 %2456
% rand('twister',s)
% si = 5;
% A = rand(si)*1.8;
% maxl = max(abs(eig(A)));
% S = 0.3*A/norm(A);

isdiscrete = 0;
verbose = 0;
compute = 1;
%S = [];
e_lon =0;
roc  = 0;
savf = 'Article';
path = '.';
namefig = 'default'


