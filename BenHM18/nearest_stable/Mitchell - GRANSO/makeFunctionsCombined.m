function [combined,N] = makeFunctionsCombined(A)

    % INITIALIZATION
    
    % First, let's only load the data once
    % load A won't work here because, when using nested fns, you can't
    % dynamically add variables to the scope.  
    % So we call a function to load the variable in its scope and assign 
    % the value to our "declared" variable.  
    A           = getA();
    is_complex  = ~isreal(A);

    % No need to compute these constants every iteration since they are...
    % fixed values.  ;-)  So let's compute them on initialization and share
    % them via nested functions and shared variables.
    if is_complex  
        n = length(A);
        N = 2*n^2;
    else
        n = length(A);
        N = n^2;
    end
     
    % RETURN the initialized function to the caller
    combined    = @computeAll;
    
    function [obj,obj_grad,ineq,ineq_grad,eq,eq_grad] = computeAll(x)
        
        % no equality constraint so these must be returned as empty 
        % we'll do here so we don't forget to later!
        eq      = [];
        eq_grad = [];
        
        % PREPARE X for both the objective and constraint
        if is_complex
            X = reshape(x(1:N/2),n,n) + 1i*reshape(x(N/2+1:N),n,n);
        else
            X = reshape(x,n,n);
        end
        
        % COMPUTE THE OBJECTIVE
        
        % f = 0.5*norm(X-A,'fro')^2; % using squared norm is not so effective
        % g = reshape(X-A,n*n,1);
        obj = norm(X-A,'fro');
        obj_G = (X-A)/obj;
        
        if is_complex == 0
            obj_grad = real(reshape(obj_G,n^2,1)); % discard imaginary rounding error
        else
            % some thought is needed to justify this.  The penalty function
            % has two terms, both real functions of complex variables: the
            % distance term and the abscissa term. Let phi denote either
            % term and let Grad denote its gradient in complex matrix
            % space.  Then we have
            %    phi(X + dX) - phi(X) \approx <Grad,dX> where <,> is the real trace
            % inner product on complex matrix space, so
            %    phi(X + dX) - phi(X) \approx Re trace Grad^* dX
            %  = (vec Re Grad)^T (vec Re dX) + (vec -i*Im Grad)^T (vec i*Im dX)
            %  = (vec Re Grad)^T (vec Re dX) + (vec Im Grad)^T (vec Im dX)
            % so setting g to the real and imaginary parts of G is the right thing.
            obj_grad = [reshape(real(obj_G),n^2,1); reshape(imag(obj_G),n^2,1)];
        end
        
        % NOW LET"S COMPUTE THE INEQUALITY CONSTRAINT 
        % X is already prepared
                
        [V,Lambda] = eig(X);
        lambda = diag(Lambda);
        % do not worry about breaking ties -- see papers for justification
        [ineq,indx] = max(real(lambda));
        v = V(:,indx); % corresp right evector (column)
        W = inv(V); % matrix of left eigenvectors
        % pert = 1e-16;
        % while isnaninf(W)
        %     % perturb V until it is numerically nonsingular
        %     V = V + pert*randn(n,n);
        %     W = inv(V);
        %     pert = pert*10
        % end
        
        w = W(indx,:); % corresp left evector
        % gradient of abscissa (ineq): see Burke and Overton, 2001 for details
        ineq_G = w'*v'; % outer product of conjugate transposes of left and right eigenvector
        
        if is_complex == 0
            ineq_grad = real(reshape(ineq_G,n^2,1)); % discard imaginary rounding error
        else
            % some thought is needed to justify this.  The penalty function
            % has two terms, both real functions of complex variables: the
            % distance term and the abscissa term. Let phi denote either
            % term and let Grad denote its gradient in complex matrix
            % space.  Then we have
            %    phi(X + dX) - phi(X) \approx <Grad,dX> where <,> is the real trace
            % inner product on complex matrix space, so
            %    phi(X + dX) - phi(X) \approx Re trace Grad^* dX
            %  = (vec Re Grad)^T (vec Re dX) + (vec -i*Im Grad)^T (vec i*Im dX)
            %  = (vec Re Grad)^T (vec Re dX) + (vec Im Grad)^T (vec Im dX)
            % so setting g to the real and imaginary parts of G is the right thing.
            ineq_grad = [reshape(real(ineq_G),n^2,1); reshape(imag(ineq_G),n^2,1)];
        end  
    end

end

function A = getA()
    load A;
end