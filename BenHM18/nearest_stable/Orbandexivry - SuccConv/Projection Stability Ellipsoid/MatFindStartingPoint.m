function [B0] = MatFindStartingPoint(A,e_lon,isdiscrete,S)
%[B0] = FindStartingPoint(A)returns a matrix which lies on the left of the
%imaginary axis.
if nargin == 1
    e_lon = 0;
    isdiscrete = 0;
    S = [];
end

%By changing here the values of the variables, you change the technique
%used for computing the starting point
Identity = 0;
Mirrored = 1;
SingularValues = 0;

if isempty(S)
    if ~isdiscrete

        if Identity
            alpha = -1;
            %             gamma = norm(A,'fro')%for debug use
            Delta = (-(A+A')+2*(alpha+e_lon)*eye(size(A)))/2;
            B0 = (A+Delta);
        elseif Mirrored
            %parameter
            alpha = -0.1;%"no-root's land" distance from the stability boundary

            %computing the eigenvalues of the matrix
            [V,ev] = eig(A);
            ev = diag(ev);

            %mirroring the unstable roots
            I = find(real(ev)>-eps+e_lon);
            if ~isempty(I)
                ev(I) = complex(-abs(real(ev(I)))-abs(e_lon),imag(ev(I)));
            end

            %shifting the roots that are too close to the stability
            %boundary
            J = find(real(ev)>alpha+e_lon);
            if ~isempty(J)
                ev(J) = ev(J)-abs(alpha+e_lon);
            end
            B0  = V*diag(ev)*inv(V);

        else
            error('Unknown type of starting point in CONTINUOUS-time in MatFindStartingPoint.m')
        end

    else
        if Mirrored
            %setting the norm of A to 1 (or to 1/gamma if the problem is
            %normed) so that all the eigenvalues of A are inside the unit
            %circle. As A is unstable, the norm of A must be greater than 1,
            %thus dividing it by its norm is acceptable
            gammac = norm(A,'fro')%for debug use
            %         B0 = A*(1+e_lon)/gamma;
            B0 = A*(1+e_lon*gammac)/gammac;

            %Another option would be to find the mirror points based on the unit
            %circle.
        elseif SingularValues
            %The idea is that the matrix with stable singular values will
            %make a good starting point for the algorithm.
            [U,S,V]= svd(A);
            s = diag(S);
            I = find(abs(diag(S))> (1+e_lon)^2);
            s(I) = (1+e_lon)^2*s(I)./abs(s(I));
            B0 = U*diag(s)*V';

        else
            error('Unknown type of starting point in DISCRETE-time in MatFindStartingPoint.m')
        end

    end
else
    B0 = S;
end
