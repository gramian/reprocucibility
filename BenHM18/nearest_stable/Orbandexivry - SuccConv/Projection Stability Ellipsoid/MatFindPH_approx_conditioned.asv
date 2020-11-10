function [H] = MatFindPH_approx_conditioned(Dp,Dq,B0,A,T,iT,n,lambda,isdiscrete,e_lon,gamma,verbose)
%Find the step H corresponding to the solution of equation (24) in
%Stable-draft.pdf. The method here is aimed at solving the system in
%"direct" way without using coordinate change or splitting of the system.
%
%Created : 13/01/2011
%Last modified : 17/01/2011

%importing the maximal radius of the ellipsoid
Stable_param;

%tolerance to apply for the convergence of the solution
tol = 1e-5;

%reconstructing some matrices
Q = T*Dq*T';
iDq = inv(Dq);
iQ = iT'*iDq*iT;
iP = iT'*iT;
iB0 = inv(B0);
P = T*T';

%evalcrit is the value measuring the convergence of the solution
crit = 0;
evalcrit = -1;
Hbar = zeros(n,n);
H = zeros(n,n);

k=1;

%% Continuous case
if ~isdiscrete
    
    %Diagonalisation using orthogonal transformation of P
    [V,L] = eig(P);
    [Vq,Lq] = eig(Q);
    Qh  = Vq*diag(   sqrt(diag(Lq)))*Vq';
    Qhi = Vq*diag(1./sqrt(diag(Lq)))*Vq';%Q^{-0.5}
    
    
    %computation of the coefficient terms
    E = ones(n,n);
    Beta = (E*L) ./ (L.^2*E + E*L.^2);
    M  = (B0-A)*P+P*(B0-A)';
    Mb = V'*M*V;
    
    %Finding mu
    mu = MatFindmu(Beta,Mb,ellips_radius,tol*1e-8);
    
    %reconstructing the step
    Xtemp = 2*mu+Beta.^2; %mu is scalar, beta is a matrix but matlab is ok with it
    Xtemp = Beta.^2./Xtemp;
    X = -Mb.*Xtemp;
    
    Htemp = X+Mb;
    Htemp = Htemp.*Beta;
    Htemp = V*Htemp*V';
    H     = A-B0+Htemp;
    

    
%% Discrete case
else
    %intialization
    scaling = (1+e_lon*gamma)^2/gamma^2;%very important do not change! (without any normalization of the problem then it is 1)
    Y =zeros(n,n);

    %computing an a transormation for P and B0*P*B0', U and V are no need
    [U,V,W,Cp,Cxx] = gsvd(T',(B0*T)');
    W = W*Cp'; %making Dp = I for this diagonalization
    iW = inv(W);
    iCp = inv(Cp);
    Dxx =  iCp'*Cxx'*Cxx*iCp;% W*Dxx*W' = B0*P*B0'
    Dxx = diag(Dxx);

    while evalcrit < crit %abs(evalcrit) > tol%AND lambda doesn't get too big (but this is implemented in the while loop

        %creating the right-hand side
        G = -(B0+Hbar - A );
        if norm(G,'fro')> eps
            G = G/norm(G,'fro');
        else
            %if the nearest matrix is inside the ellipsoid (the original matrix is
            %stable) go to that point
            H = Hbar;
            return
        end
        Rhs = iW*Q*G*iP*iB0*Q*iW';

        %solving O(n^2) small 2*2 subsystems (upper triangular part of the
        %matrix)
        %diagonal terms
        for i = 1:n
            if abs(Dxx(i)) ~= scaling
                Y(i,i) = complex(real(Rhs(i,i))/(scaling+Dxx(i)),imag(Rhs(i,i))/(scaling-Dxx(i)));
            else
                %A bug might occur since the next equation does not detect an
                %impossible equation if the complex or real part of Rhs is not
                %zero when Dxx(i) =scaling or -scaling respectively
                Y(i,i) = (scaling+Dxx(i))/2*complex(real(Rhs(i,i)),       0      )+...
                    (scaling-Dxx(i))/2*complex(       0      ,imag(Rhs(i,i)));
            end
        end

        %off-diagonal terms
        for j = 2 : n
            for i = 1 : j-1

                %building the 2*2 subsystem
                Mre = [scaling Dxx(i); Dxx(j) scaling];
                Nre = real([Rhs(i,j); Rhs(j,i)]);
                Mim = [scaling -Dxx(i); -Dxx(j) scaling];
                Nim = imag([Rhs(i,j); Rhs(j,i)]);

                %solving it
                yre = Mre\Nre;
                yim = Mim\Nim;

                %saving in a matrix at position (i,j) and (j,i)
                Y([i+n*(j-1) j+n*(i-1)]) = complex(yre,yim);
            end
        end


        %reconstructing the solution
        Hprev = H;
        H = W*Y*iW*B0;

        %evaluating lambda
        YDxx = Y*diag(Dxx);
        BH = W*(YDxx+YDxx')*W';
        Fp1 = real(trace(iQ*BH*iQ*BH'));
        WY = W*Y;
        Fp2 = real(trace(iQ'*2*WY*diag(Dxx)*WY'));
        FpXHH = Fp1+Fp2;
        lambda = sqrt(4*ellips_radius)/sqrt(FpXHH);

        %setting the length of H so that it is on the boundary of the
        %ellipsoid
        H = lambda/2*H;

        if any(eig((B0+H)*P*(B0+H)'-P)>0)
            'coucou'
        end
        %tstar is the oefficient for an optimal step, it can be positive or
        %negative, and bigger or smaller than one.
        normHHbar = norm(H-Hbar,'fro');
        tstar = -real(trace((B0+Hbar-A)*(H-Hbar)'))/normHHbar^2;

        relax = 0;
        if tstar > 1
            tstar = 1;
        elseif tstar < 0
            tstar = 0;
        elseif tstar < tol*normHHbar
            relax = 1;
        end

        %         %DEBUG : This code is to know if tstar gets B0+Hbar-A+tstar*(H-Hbar)
        %         %        to its minimum
        %         tint = [tstar-5:0.2:tstar+5];
        %         for i = 1 : length(tint)
        %             yint(i) = norm(B0+Hbar+tint(i)*(H-Hbar)-A,'fro');
        %         end
        %         hold on;
        %         plot(tint,yint,'b',tstar,norm(B0+Hbar+tstar*(H-Hbar)-A,'fro'),'or')
        
        %The point B0+Hbar is the point where we consider the derivative to
        %the distance to A. It is allowed to get outside the ellipsoid.
        Hbarprev = Hbar;
        Hbar = Hbar + tstar*(H-Hbar);

        %when the length of the radius of the ellipsoid doesn't change much
        eprev = evalcrit;
        evalcrit =norm(H-Hprev,'fro')/norm(H,'fro');
        evalcrit = (norm(B0+Hbarprev-A,'fro')-norm(B0+Hbar-A,'fro'))/norm(Hbar,'fro');
        evalcrit =norm(Hbar-Hbarprev,'fro')/norm(Hbar,'fro');   
        evalcrit =(norm(Hbar,'fro')-norm(Hbarprev,'fro'))/norm(Hbar,'fro');  
        evalcrit = real(trace((B0+Hbarprev-A)*(bpram*H-Hbarprev)'));
        
        
        if relax == 1
            crit = crit -tol*normHHbar/10;
        end
%         if abs(eprev - evalcrit) < eps
%             'debug'
%         end
        
        k= k + 1;
        %         if norm(H-Hbar,'fro')/norm(H,'fro') < tol
        %             H2 = H;
        %         end
        %         if norm(H-Hprev,'fro')/norm(H,'fro') < tol
        %             H3 = H;
        %         end
        %         dt = norm(B0+Hbar-A,'fro')
        %         dl = norm(B0+H-A,'fro')
        %
        %         evalcrit1 = norm(tstar*(H-Hbar),'fro')
    end
    k
    crit
    XPH = B0*P*Hbar';
    BH = XPH+XPH';
    Fp1 = real(trace(iQ*BH*iQ*BH'));
    Fp2 = real(trace(iQ*2*Hbar*P*Hbar'));
    FpXHH = Fp1+Fp2;
    lab = sqrt(4*ellips_radius)/sqrt(FpXHH);
        
    if norm(B0+Hbar-A,'fro') < norm(B0+lab/2*Hbar-A,'fro')
        H = Hbar;
    else
        H = lab/2*Hbar;
    end
%%
end
%WARNING: This routine was built for a Dp = I ONLY
