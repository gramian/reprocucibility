function [b0] = PolySPReflection(A,e_lon,isdiscrete,coeff)
%[b0] = PolySPReflection(A,e_lon,isdiscrete) creates a stable polynomial b0
%based on reflection coefficients. For a discrete polynomial, reflection
%coefficients will all be less than one if the polynomial is stable. This
%method thus compute reflection coefficients of a = poly(A), and force them
%below 1 whenever they are greater.
%
%In case, A represents a continuous time polynomial, a Moebius transform to
%the discrete setting is achieved so that the problem is solved on an
%equivalent discrete polynomial and then brought back to a continuous time
%setting.
%
%The number of roots out of the unity circle is given by the number of
%negative 'Pi'. Pn is given by, for 1<=i<=n-1 (where n is the degree),
% Pi = prod(1-abs(k(j))) for j = [i : n-1]; 
%
%Input : A         : a polynomial/companion matrix representing a polynomial.
%        e_lon     : a scalar that gives a distance to stability so that b0
%                    is not too close to the boundary (of the unity circle)
%        isdiscrete: 1 if the problem is discrete, 0 otherwise.
%
%History : 04/01/2011 : Creation
%          19/01/2011 : StepDown and StepUp procedure modified for complex
%                       coefficients polynomials
%          20/01/2011 : StepDown and StepUp return a line instead of a
%                       column
%          25/01/2011 : Changing the instable coefficients, k(I), are
%                       modified to a stable state.
%          30/09/2011 : A can be a polynomial or a companion matrix
%Last modified : 30/09/2011 

%% Preprocessing

%parameter for the method : coeff enables you to control how the instable
%reflection coefficients will be modified. Typically, the reflection
%coefficients are mirrored w/r to the unit circle, and multiplied by coeff
%so that it is not exactly mirrored
% coeff = 0.99;

if nargin < 4;
    %coeff is best when nonpositive (for the distance to 'a') and close to
    %zero (for staying away from the instable region).
    coeff = 0;
end
if nargin < 3;
    isdiscrete = 1;
end

if nargin < 2;
    e_lon = 0;
end

%vector of coefficients of the polynomial
[nc nl] = size(A);
if nc==nl
    a = poly(A);
elseif nc == 1;
    a = A.';
elseif nl ==1;
    a = A;
else
    error('Error in PolySPReflection : A should be a square companion matrix or a polynomial');
end

if ~isdiscrete
   a = PolyBijectionDC(a); 
   a = a./a(1);
end


%% Processing
%The vector 'a' is the vector of coefficients of the polynomial A_M(z)
%which has the expression :
%                   A_M(z) = \sum_{i=0}^M a_i z^{-i)
%Note that the roots of A_M are identical to the roots of 'A'.
%Given that indices start with 1 in Matlab, the sum is shifted of one.
%In order to find the reflection coefficients, we implements the stepdown
%procedure of the book of Markel and Gray. 
%We then correct the reflection coefficients so that |k_m| < 1 for all m.
%Finally, and if the reflection coefficients have been modified we
%reconstruct a stable polynomial from the new reflection coefficients.

%k is a vector containing the reflection coefficients of the polynomial
k = StepDown(a);

%modification of the reflection coefficients so that the modulus of the co-
%mplex reflection coefficient is now 1/kmodul, but keeps the same structure
kmodul = abs(k);
I = find(kmodul >= (1-eps));
k(I) = (1- coeff)*k(I)./kmodul(I);% for tests : ./kmodul(I).^2; 


%Applying the stepup procedure to obtain a stabilized polynomial from the
%modified reflection coefficients
if ~isempty(I)
    b0 = StepUp(k);
else
    b0 = a;
end

%achieving the inverese mapping if we are in the continuous time setting
if ~isdiscrete
   b0 = PolyBijectionDC(b0); 
end

%the polynomial must become monic
b0 = b0./b0(1);

%ordering by ascending power instead of descending, 
%b0 should be a column
b0 = b0(end:-1:2).';
end

%%
function [k] = StepDown(a)
%Computes the reflection coefficients of a polynomial;
%Ref : "Linear Prediction of speech", J.Markel and A.Gray,1976,Springer-Verlag
%      "Immitance-type three-term schur and levinson recursions for quasi-
%      toeplitz complex hermitian matrices", Y.Bistritz, H. Lev-Ari and T.Kailath, 
%      SIAM J.Matrix Analysis and applications,Vol 12, pp 497-520

M = length(a)-1; % M is thus the degree of a
k = zeros(M,1);
u = a;

for m = M:-1:1
    %initializing for the current iteration
    v = zeros(1,m);
    km = u(m+1);
    k(m) = km;
    den = 1-abs(km)^2;%Pay attention to the fact that km might be 1.
    
    %step down procedure
    for i = 1:m  
        v(i) = (u(i)-km*conj(u(m+2-i)))/den;
    end
    u = v;
end

end


%%
function [b] = StepUp(k)
% this function is the inverse function of StepDown
M = length(k);%degree of the polynomial
v = 1;%base polynomial

for m = 1:M
    km = k(m);
    u = zeros(1,m+1);
    for i = 1:m+1
        if i == 1
            u(i) = v(i);
        elseif i == m+1
            u(i) = km;
        else
            u(i) = v(i)+km*conj(v(m+2-i));
        end
    end
    v =u;
end
b = u;
end

    



