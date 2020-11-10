function [H] = MatFindPH_approx_updated(Dp,Dq,B0,A,T,iT,n,lambda,isdiscrete,e_lon,gamma,verbose)
%Find the step H corresponding to the solution of equation (24) in
%Stable-draft.pdf. The method here is aimed at solving the system in
%"direct" way without using coordinate change or splitting of the system.

%WARNING: This routine has been built for a Dp = I ONLY
%this warning is helful if find_PQ has been changed inadvertently with a
%Dp ~= 1
% if any(Dp ~= speye(n))
%     error('FindPH_approx:DpNotIdentity','Dp is not identity and FindPH_approx should therefore not be used')
% end

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
la = 1;
%%
if ~isdiscrete
    
%Precomputation  
    %computing an orthogonal decomposition for P(ONLY accurate for Dp = I)
    [U,S,V] = svd(T);

    %P = T*T' if Dp = I, and P^2 = (T*T')*(T*T')
    S2 = S*S';
    S4 = S2*S2;

%Optimizing on the kernel    
    %computing the optimal dual variable solution of the optimization
    %problem on the kernel of the Lyapunov operator
    %         Y = lyap(2*Ps,-((A-B0)*P+P*(A-B0)'));

    %We first construct an equivalent to M = (A-B0)*P+P*(A-B0)';
    Num = U'*(A-B0)*U*S2;
    Den = ones(n,n)*S4;%will never contains zeros
    Den = Den+Den.';
    Mu  = (Num+ Num');
    Y = U*(Mu./(2*Den))*U';
    Yc = T'*Y*T;
    YP = Y*(T*T');

%Solving the two lyapunov equations    
    %iterating to find the right H
    while (evalcrit < crit) %the tol should vary with the size of the ellipsoid
        
        %it now remains to solve the problem of equalling the two
        %derivatives
        G = (-2*Yc + T'*Hbar*iT')/la; % the division by lambda could be removed but it enhances the numerical stability
        if norm(G,'fro')> eps
            G = G/norm(G,'fro');%it doesn't change anything to norm it, the derivative is still the same
        else
            %if the nearest matrix is inside the ellipsoid (the original matrix is
            %stable) go to that point.
            H = Hbar;
            return
        end    
        Ht = -1/(2)*Dq*G*Dq ;%the transformation needed for hbar is different that the one for h
        H = T*Ht*iT;

     %computing the correction to apply in order to have (H+Dta)*iP
     %hermitian
        %the skew-hermitian part of Dta
        Hr = T'*H*iT';
        Dta_a = -0.5*(Hr-Hr');

        %Finding the hermitian part of Dta
        Num = V'*Dta_a*V*S4;
        Mv = Num + Num';
        Dta_h = -V*(Mv./Den)*V';

        %Reconstructing Dta, now (H+Dta)iP should be hermitian
        Dta = iT'*(Dta_a+Dta_h)*T';

        %adding the correctino to H
        H = H + Dta;
        
     %Computing lambda   
        %computing the scaling factor so that H is on the boundary
        HppH = H*P+P*H';
        FpXHH = real(trace(iQ*HppH*iQ*HppH'));
        la = sqrt(4*ellips_radius)/sqrt(FpXHH);
        
        %adjusting H
        H = la/2*H;
        
     %Optimizing on the line (X+H, X+Hba)   
        %tstar is the oefficient for an optimal step, it can be positive or
        %negative, and bigger or smaller than one.
        normHHbar = norm(H-Hbar,'fro');
        tstar = -real(trace((-2*YP+Hbar)*(H-Hbar)'))/norm(H-Hbar,'fro')^2;

        relax = 0;
        if tstar > 1
            tstar = 1;
        elseif tstar < 0
            tstar = 0;
        elseif tstar < tol*normHHbar
            relax = 1;
        end

      %Saving the step and the previous one  
        Hbarprev = Hbar;
        Hbar = Hbar + tstar*(H-Hbar);
%         d = norm(B0+Hbar-A,'fro');
        Hprev = H;
        
      %Evaluating the criterion, the derivative at the point X+Hbarprev
      %should not intersect with an interior ellipsoid which is "bpram"
      %percent of the original one. (As it requires to compute H_{k+1} we
      %check for the previous step)
        %computing the norm of the change in step
        %         evalcrit = abs((norm(Hbar,'fro')-norm(Hbarprev,'fro')))/norm(Hbar,'fro');
        %         evalcrit = abs(norm(B0+Hbarprev-A,'fro')-norm(B0+Hbar-A,'fro'))/norm(B0+Hbar-A,'fro');
        evalcrit = real(trace((-2*YP+Hbarprev)*(bpram*H-Hbarprev)'));
        
        
        if relax == 1
            crit = crit -tol*normHHbar/10;
        end

        %         evalcrit = norm(H-Hprev,'fro')/norm(H,'fro');
        k = k +1;
   
    end

    HppH = Hbar*P+P*Hbar';
    FpXHHb =real(trace(iQ*HppH*iQ*HppH'));
    lab = sqrt(4*ellips_radius)/sqrt(FpXHHb);

    if norm(B0+Hbar-A,'fro') < norm(B0+lab/2*Hbar-A,'fro')
        H = Hbar;
    else
        H = lab/2*Hbar;
    end
    k;
    H = (A-B0-2*YP)+H;%Xh + Hbar
    %     res = 2*real(trace(iQ*(H*P+P*H)*iQ*P*H'))
%%
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
