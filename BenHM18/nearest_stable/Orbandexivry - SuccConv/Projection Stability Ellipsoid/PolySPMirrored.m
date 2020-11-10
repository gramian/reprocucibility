function [b0] = PolySPMirrored(A,e_lon,isdiscrete)
%[b0] = PolySPReflection(A,e_lon,isdiscrete) creates a stable polynomial b0
%based on the reflection of the unstable roots of the polynomial. 
%
%The computation of the roots is achieved using the matlab method "root".
%Unstable roots are then mirrored with respect to the stability boundary
%(the unit circle in discrete time, the imaginary axis in continuous time).
%
%
%Input : A         : a companion matrix representing a polynomial or a
%                    polynomial
%        e_lon     : a scalar that gives a distance to stability so that b0
%                    is not too close to the boundary
%        isdiscrete: 1 if the problem is discrete, 0 otherwise.
%
%History : 01/07/2011 : Creation
%          28/03/2012 : Modifies the file so that it can deals with
%                       matrix and polynomial inputs (needed for FSP)
%Last modified : 28/03/2012

%% Preprocessing

%parameter : default_dist specifies the distance at which the roots that
%are marginaly stable should be put
default_dist = 0.1;

[s1,s2] = size(A);

if s2 == 1
    A = A.';%if A is a polynomial it must be a line (on the contrary of FSP)
end
if s1 == s2
    %vector of coefficients of the polynomial
    a = poly(A);
elseif any([s1 s2] ==1)
    a = A;
else
    error('PolySPMirrored:wronginputsize','The input is not a vector or a square matrix')
end


%% Processing

%computing the roots of the polynomial in O(n^3) operations
allroots = roots(a);

if ~isdiscrete
    
    %finding the roots that are unstable
    I = find(real(allroots) > 0+e_lon);
    
    %finding the roots that are on the stability boundary
    J = find(abs(real(allroots))-e_lon < default_dist);
    
    %mirroring them with respect to the imaginary axis
    allroots(I) = e_lon - (conj(allroots(I))-e_lon); 
    allroots(J) = allroots(J)- default_dist;
    
else    
    rootsmodul = abs(allroots);
    %finding the roots that are unstable
    I = find(rootsmodul > 1+e_lon);
    
    %finding the roots that are on the stability boundary
    J = find(rootsmodul-1-e_lon < default_dist);
    
    %mirroring (in the sense inversely proportional) them with respect to
    %the imaginary axis
    allroots(I) = allroots(I).*(1+e_lon).^2./(rootsmodul(I).^2); 
    allroots(J) = allroots(J)*(1-default_dist);
    
end


%reconstructing a stable polynomial
b0 = poly(allroots);
b0 = fliplr(b0(2:end)).';
