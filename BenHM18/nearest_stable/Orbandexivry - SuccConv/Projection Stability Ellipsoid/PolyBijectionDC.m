function [q] = PolyBijectionDC(p)
%[q] = PolyBijectionDC(p) makes a change of coordinate so that the
%left half-plane is in bijection with the unit disk. The change of variable
%is its auto inverse.
%Input : p : A vector containing the n+1 coefficients of a polynomial of
%            order n in descending order.
%Output: q : The bijection of p in the unit disk, or its inverse if p 
%            already lies in the unit circle 
%Created       : 03/12/10
%Last modified : 06/12/10

% if nargin == 0
%     p = [-2 3+j 1 4-5j 6 ];
% end

n = length(p)-1; %degree of the polynomial

%computing the powers of (z-1) and (z+1) (only those are needed)
lmid = floor((n+1)/2);%lower middle element
umid = ceil((n+1)/2); % upper middle element
gap = umid-lmid;

%allocating space for the table of the power of (z-1) until the row
%[ceil((n+1)/2)-1] and the power of (z+1) from floor((n+1)/2) afterwards
%until (z+1)^n. It has n+1 or n+2 lines depending on odd or even n. If
%there is n+2 lines the "middle" lines (at umid and lmid+1) have
%polynomials of equal powers.
powofz = zeros(n+1 + gap,n+1);
powofz(1,n+1) = 1; %unit polynomial :(z-1)^0

%base polynomials
polm = [1 -1]; %(z-1)
polp = [1  1]; %(z+1)

%initiation
apol = 1;%unit polynomial
lpol = 1;%length of apol
for i = 1 : umid-1 % the -1 is because there is already one power inside the table
    if gap
        preapol = apol;
    end
    lpol = 2+lpol-1;
    apol = conv(apol,polm);
    powofz(i+1,n+2-lpol:n+1) = apol;
end

%transforming to get the powers of (z+1)
if gap
    %in case n is even, we have to create two polynomials of degree lmid
    apol = preapol;
    lpol = lpol-1;
end
apol = abs(apol);

for i = lmid : n%lmid is always the first power of (z+1) to be computed whenever n is odd or even.
    lpol = 2 + lpol-1;
    apol = conv(apol, polp);
    powofz(i+1+gap,n+2-lpol:n+1) = apol;
end

%computing the polynomials (z+1)^(n-i)*(z-1)^(i)
prod = zeros(n+1,n+1);
for i = 0: umid-1

    
    %We start from the middle coefficient and we go to the right, in the
    %table of powers, it means starting at the middle and goind to the
    %extremes.
    %taking the two polynomials from the table
    pol1 = powofz(umid-i,n+2-umid+i:n+1);%taking a polynomial of length j (z-1)^{mid-j-1)
    pol2 = powofz(umid+1+i,n+2-(umid+1+i-gap):n+1);
    
    %multiplying the two polynomials 
    poltemp = conv(pol1,pol2);
    
    %saving their values 
    prod(umid-i,:) = poltemp;
    prod(umid+(1-gap)+i,:) = poltemp;
end

%mirroring the coefficients with their signs for computing the rest of the
%tabular (z-1)^(n-i)*(z+1)^(i) (the coefficients z+1 and z-1 have exchange
%positions). For p is ordered in ascending order, each line must correspond
%to it hence (z+1)^n will be at the beginning of the tabular and (z-1)^n
%will be at the end.
J = [2:2:n+1];
prod(umid+1:n+1, J) = -prod(umid+1:n+1, J);

%multiplying with the coefficient of p and summing to obtain q
q = p*prod;


    
    
    