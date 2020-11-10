%discrete
X = Bi;
iT = invT;
P = T*Dp*T';
iP = iT'*inv(Dp)*iT;
iQ = iT'*inv(Dq)*iT;

XiQ = X*iQ;
PXiQ = P*XiQ;
PAP = P+PXiQ*X'*P;

2*real(trace((PXiQ*H'*PXiQ + PAP*H*iQ)*H'))

iY = [iP+X*iQ*X' -X*iQ; -iQ*X' iQ];
R = [zeros(size(P)) P*H; H'*P zeros(size(P))];

real(trace(iY*R*iY*R'))


%continu
bP= 2*iQ*(H*P+P*H')*iQ*P;

real(trace(bP*H'))