function [b0] = PolySPHenrici(A,e_lon,isdiscrete)
%[b0] = PolySPHenrici(A,e_lon,isdiscrete) creates a stable polynomial b0
%based on Henrici's method. For a discrete polynomial, Henrici's method
%computes the fast fourier transform of the polynomial and links it to the
%laurent serie approximation of the same polynomial. This link enables an
%incomplete factorization of the polynomial in a purely stable polynomial
%and a purely unstable polynomial. Mirroring the coefficients of the
%unstable part gives us a final stable polynomial.
%
%In case, A represents a continuous time polynomial, a mapping to the
%discrete setting is achieved so that the problem is solved on an
%equivalent discrete polynomial and then brought back to a continuous time
%setting.
%
%Input : A         : a companion matrix representing a polynomial.
%        e_lon     : a scalar that gives a distance to stability so that b0
%                    is not too close to the boundary (of the unity circle)
%        isdiscrete: 1 if the problem is discrete, 0 otherwise.
%
%History : 05/12/2011 : Creation
%          05/01/2011 : Making of the header
%Last modified : 05/01/2011 


%tolerance for the Fourier approximation
tol = 1e-14;

%vector of coefficients of the polynomial
a = poly(A);
if ~isdiscrete
   a = PolyBijectionDC(a); 
   a = a./a(1);
end

%line vector with the coefficients
%     a = [ -A(:,end);1].';
%     a = a(lena:-1:1);

%length of a
lena = length(a);

%computing the approximation of the laurent coefficient of f(z) and the order
%of the fourier series so as to have acceptable error
[n,am,k] = PolyStart_LaurentCoef(a,tol);

%now that the sum of the powers of the roots are known. we can try to
%compute the coefficients of r(z) = prod(inside roots). This is done
%using multiple imbricated newton scheme.

%number of sums of power s^1 s^npow
npow = length(am)-2;
m    = [1:npow];

%constructing the sums of power
spow = am(3:end)./m;%am(2) correspond to the sum of the roots to the power 1
lenq = npow+1;
Q    = -[spow(npow:-1:1) 0];%Q(end) is the indepedent term which is 0 for a non unit.

%newton algorithm to get qi(z) whose has the same coefficients as pi(z)
%but the inverse roots
nite = ceil(log2(k))+1;% the "+1" is added so that we can check that the right solution is obtained
W    = 0;

for h = 1:nite
    exactnum = 2^h;%the number of accurate decimals at the end if the current iterate
    lenw   = exactnum/2;
    id     = [1: exactnum];

    %constructing the first term
    term1  = W ;
    term1(lenw) = term1(lenw) + 1;% = [zeros 1]+W

    %constructing the second term
    LogW   = PolyStart_Log(W);
    index  = lenq-exactnum +id; % npow will always be greater than exactnum since k is small
    term2  = LogW - Q(index);
    stepup = -conv(term1,term2);

    %truncating to keep only the exactnum last terms
    index  = lenw-1 + id; %stepup has a size exactnum + lenw -1 (=3*lenw -1)
    stepup = stepup(index);

    %adding to W (minus sign already in stepup)
    index = lenw + [1:lenw];
    stepup(index) = stepup(index)+ W;

    %allocating the result
    W = stepup;
end

%checking the quality of the solution
if abs(W(2*lenw-k-1))> tol
    warning('The localisation  of the stable roots of the polynomial might lack precision,\nThe user cannot address this warning.')
end

%constructing q1 that has roots inverse to those of p1
index = 2*lenw-k-1 +[1:k+1];
q1    = W(index);

%Moreover in a Companion matrix we consider the polynomial to be monic
%(as is q1 after the next commented line) and we do not write the 1 in
%the companion matrix. To get the right roots of q1 one should use the
%following line :
q1(k+1) = q1(k+1) + 1;
%     %but we delete this index so that we have
%     q1 = q1(1:k);

%     %we have to return q1 since the matlab convention and ours are
%     %opposite, we consider coefficients in ascending order and matlab
%     %consider it in descending order.
%     b0= q1.';

%For information only : the last step to get a polynomial with the
%right roots would be :
%Constructing p1 which has the same coefficients as q1 but in reverse
%order the roots of p1 contains the roots inside the circle centered at
%0 and of radius 1
%
p1 = q1(k+1:-1:1);

%building the inverse of p1 at a precision equivalent to the length of
%p2
lenp2 = lena-k;%=lena-(k+1)+1, lenp2 can never become zero since k counts the number of roots and lena compute the number of coefficients
p1inv = PolyStart_inv(p1,lenp2);
p2 = conv(p1inv,a);

%truncating p2 to the right length so that lenp2+lenp1-1 = lena
index = lena-1 +[1:lenp2];%=lenp2+lena-1 - lenp2 +[1:lenp2]
p2 = p2(index);

%constructing a stable polynomial
%the comments reflects the fact that we are only dealing with discrete
%problem right now. 
% if isdiscrete
    %in the discrete case you take a polynomial whose roots are the
    %mirrored roots of p2 to get p2 stable (and hence b0).The mirrored
    %roots when the boundary is the unit disk are the conjugate inverse
    %roots.
    p2 = conj(p2(lenp2:-1:1));
% elseif lenp2 > 1
%     %if we are in the continuous case if p2 has at least one root, we
%     %take the opposite
%     p2  = -p2;%taking the inverse does not suffice, you have to think
%     which signs should change.
% end

%building the stable polynomial that will be used as a starting point
sol = conv(p1,p2);

%if the problem is continuous, applying the transformation a second time
%gives us the approximation of the original polynomial. Sol must be is
%descending order
if ~isdiscrete
    sol = PolyBijectionDC(sol);
end

%the polynomial must become monic
b0 = sol./sol(1);

%ordering by ascending power instead of descending
b0 = b0(end:-1:2).';

end



%%
function [n,am,insidroot] = PolyStart_LaurentCoef(p,tol)
% [n,am, insidroot] = PolyStart_LaurentCoef(p,tol) computes the approx-
% imation of the Laurent coefficient of the function f(z) = p'(z)/p(z).
% The algorithm presented here is based on the section 3.1 of the paper
% from Henrici, "Fast fourier methods in computational complex analysis",
% SIAM 1979.
% The method iterates on n, until am(2) is within distance 'tol' of an
% integer.
% Input : p       : the coefficient of the polynomial p ordered by
%                   decreasing power
%         tol     : The tolerance for the exactness of the approximation,
%                   it is measured on the a_{-1} which should be an integer
%                   that refers to the number of zeros inside the circle
%                   centered on 0 and of radius 1.
% Ouput: n        : The number of roots of unity used to compute the fourier
%                   coefficients used as approximants. The smallest the
%                   tolerance, the highest n. Initial value : n = 20.
%        am       : a vector of length n containing the approximation of the
%                   Laurent coefficients, am(1) = am_{0}, am(2) = am_{-1},
%                   am(n) = am_{-(n-1)}
%        insidroot: the number of roots inside the complex circle centered
%                   at 0 and of radius 1.
%
% Created : 4/11/2010
% Last modification : 4/11/2010

%derivative of p
pder = polyder(p);

%degree of the polynomial
deg = length(p);

%initial n
n = 20;
closeness = 1;

%the while measure the closeness of a_{-1} to an integer, a_{-1}
%indicates the number of zeros inside the circle
while closeness > tol

    w = exp(i*2*pi/n);%complex

    k = [1:n];
    %constructing f(w^k)
    num = polyval(pder,w.^(k-1));
    den = polyval(p   ,w.^(k-1));
    f = num./den;

    %constructing the approximation of the first laurent coefficients a_{-1}
    am = ifft(f);

    %the second element of am corresponds to a_{-1}
    insidroot = round(am(2));
    closeness = abs(insidroot-am(2));

    n = n + 10;
    if n >500
        n = n + 990;
    end
    if n > 7e4
        closeness= 0;
        warning('At least one of the roots has a modulus close to 1. The factoring might lack accuracy.')
    end
end
end

%%
function [LW] =  PolyStart_Log(W)
% [LW] =  PolyStart_Log(W) computes the logarithm of a given polynomial
% W of order N. The logarithm is computed as the integral of a its
% derivative. Given that the derivative is a rational function of two
% polynomials, a Newton scheme is applied in order to solve the
% fraction.
% Input : W  : A non-unit polynomial (0 is a root of the polynomial)
%              ordered by decreasing power.
% Output: LW : LW the polynomial approximating the logarithm of W and
%              whose order is twice the order of W.
%
% Created : 4/11/2010
% Last modification : 5/11/2010

lenw = length(W);

%derivative of W
Wder = polyder(W);

%constructing 1+W
Wone = W;
Wone(lenw) = Wone(lenw)+ 1;%=W +[zeros 1] since W is a non-unit

%computing the inverse of 1+W, Y has length 2*lenw
[Y] = PolyStart_inv(Wone,2*lenw);

%computing the derivative of the log
LWder = conv(Wder,Y);

%truncating the development to the desired length
index = [1:2*lenw-1];%number of terms that we keep in the polynomial approximations
indexcut = (lenw-1) + index ;%=(lenw-1) +(2*lenw) -1 -(2*lenw-1) +index
LWder = LWder(indexcut); %we keep only 2*l-1 terms given we wan to integrate to get 2*l terms

%coefficient for the final integration and integration
LW = polyint(LWder,0);
end


%%
function [Y] = PolyStart_inv(P,len)
% [Y] = PolyStart_inv(P) inverts the polynomial P by executing a newton
% scheme on P. The ouput Y is a polynomial whose length is specified by
% len, except if P is a scalar. The inverse is computed in the sense that
% I = conv(Y,P) = 1. The result I is only accurate to the order len-1 (or
% 'len' coefficients)
%
% Input : P   : A vector of the coefficients of a polynomial ordered by
%               decreasing power. P must be a unit, i.e P(end) ~=0.
%         len : The desired length (precision) of Y.
% Output: Y   : The polynomial inverse of P
%
% Created : 05/11/2010
% Last modification : 05/11/2010

lenp = length(P);
if P(lenp) == 0
    error('P must be a unit to be inversed, i.e. P(end) ~=0, error in PolyFindStartingPoint.m')
end

if lenp == 1
    Y = 1/P;
    return
end

%initialisation of the newton scheme
nite = ceil(log2(len)); %number of iteration that are needed to reach the desired accuracy
Y = 1/P(lenp);%initial value, Y is the inverse of P
lp2 = 1;%length of Y at his stage, it initialized because nite can be zero if len is 1 in which case lp2 has to exist

for h = 1 : nite
    %length of Y at the end of iteration
    exactnum = 2^h;
    id = [1 : exactnum];

    %first convolution
    lp1   = lenp+exactnum/2-1;%=length(prod1);
    prod1 = conv(P,Y);
    prod1(lp1) = prod1(lp1)-2;

    %truncating the first product so to keep the exact coefficients
    %only
    if lp1 > exactnum
        index = lp1-exactnum+id;
        lp1 = exactnum; %the length of prod1 is modified (shorter)
    else
        index = [1:lp1];
    end
    term1 = -prod1(index);

    %second product and truncation
    lp2   = exactnum/2 + lp1 -1;
    prod2 = conv(Y,term1);

    %Truncating so to get Y with accurates numbers only
    if lp2 > exactnum
        index = lp2-exactnum +id;
        lp2 = exactnum;
    else
        index = [1:lp2];
    end
    Y = prod2(index);

end

%truncating the final iterate to a length 2*lenp if it is too long
if lp2 > len
    index = lp2-len + [1:len];
else
    index = [1:lp2];
end
Y = Y(index);

end
