function [Dp,Dq,T,invT] = MatFindPfromQrandom(B,S,n,e_lon,isdiscrete,verbose)

B = full(B);
if isempty(S)
    S = 100*rand(n);
end

if ~isdiscrete
    %Computing L_q^{-1} solution of the lyapunov equation
    %-nI = B_i^*Q^{-1}+Q^{-1}B_i
    Y = B-eye(n)*e_lon;

    %computing the cholesky factor Lq
    Lq = chol(S*S');

    %determining Lp the cholesky factor of P
    Up = lyapchol(Y,S);

    if verbose
%         P = Up'*Up
%         Q = Lq'*Lq

%         eigP = eig(P)
%         eigQ = eig(Q)
    end
%%
else

    
    %As the problem is normalized, you want the eigenvalues to lie inside a
    %circle of radius 1/gamma. Thus the following equations/change of
    %variables apply in the finding of P and Q.
    Y = B*1/(1+e_lon);

    %Q = S*S';
    Lq = chol(S*S');
    
    %determining Lp the cholesky factor of P
    Up = dlyapchol(Y,1/(1+e_lon)*S);

    if verbose
%         P = Up'*Up
%         Q = Lq'*Lq

%         eigP = eig(P)
%         eigQ = eig(Q)
    end
%%
end

%We now need to diagonalize P and Q, as we have them under a Cholesky form,
%finding thier diagonalized form is similar to solve a generalized singular
%value decomposition.
[U,V,T,Cp,Cq] = gsvd(Up,Lq);

T = T*Cp';
iCp = inv(Cp);
Dp = speye(n);
Dq = sparse(iCp'*Cq'*Cq*iCp);
% Dp = sparse(Cp'*Cp);
% Dq = sparse(Cq'*Cq);

invT = inv(T);
