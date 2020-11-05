function [A,B,C] = isp(N,IO,r)
% isp (inverse sylvester procedure)
% by Christian Himpe, 2016 ( http://gramian.de )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%*
    if(nargin==3)
        rand('seed',r);
        randn('seed',r);
    end;

    % Gramian Eigenvalues
    a = 1e-1;
    b = 1e+1;
    WX = -diag( a*((b/a).^rand(N,1)) );

    % Gramian Eigenvectors
    [Q,R] = qr(randn(N,N));

    % Input and Output
    B = randn(N,IO);
    C = B';

    % Solve System Matrix
    A = sylvester(WX,WX,B*C) - sqrt(eps)*speye(N);

    % Unbalance
    A = Q'*A*Q;
    B = Q'*B;
    C = C*Q;
end
