function prob = makeProblem(U0,X0,X1,Y1,user_opts) 
%   makeProblem:
%       Takes the problem data and constructs the necessary functions and
%       metadata for running GRANSO and recovering the solution(s) it
%       finds.  This function also provides the code for computing the
%       gradients of the objective and inequality constraint functions.
%
%   USAGE:
%       p = makeProblem(U0,X0,X1,Y1,user_opts)
%
%   REQUIRED INPUTS:
%       U0,X0,X1,Y1
%           matrices defining the optimization problem
%
%   OPTIONAL PARAMETERS: 
%       user_opts           
%           A struct containing any or all of the following fields:
%
%           .discrete_time          [ boolean | {true} ]
%               True if the system is discrete time, and thus the stability
%               measure should be the spectral radius minus 1. False if the
%               system is continuous time, in which case the stability
%               measure is then the spectral abscissa.
%
%           .frobenius_obj          [ boolean | {true} ]
%               Whether or not to use the Frobenius norm squared or the
%               spectral norm for the objective function.
%
%           .regularization         [ nonnegative real | {0} ]
%               When this is zero, no regularization is used.  When this is
%               positive, the following is added to the objective function:
%       
%                   opts.regularization * regularization_term 
%            
%               opts.regularization_type chooses the specific term added.
%
%           .regularization_type    [ positive integer | {1} ]
%               Chooses the specific regularization term to add:
%               1:  spectral norm squared of variables 
%               (new types may be added later for larger integer values)
%
%           .require_stable         [ boolean | {true} ]
%               Set up a constrained problem where a solution is infeasible
%               if F is not deemed stable (as governed by the relevant
%               options).
%   
%           .stability_margin       [ nonnegative real | {0} ]
%               Margin for stability.  By default, a solution is stable if
%               its stability constraint is less than 0 (e.g. spectral
%               abscissa or spectral radius minus 1).  If stability_margin
%               is set to a positive value, then a solution is only deemed
%               stable if its stability measure is less than
%               -stability_margin.
%
%   OUTPUT:
%       prob  
%           Struct of data and function handles for the problem.  The first
%           fields replicate the values set by user_opts (or if not
%           specified by the user, the default values).  The subsequent
%           fields are:
%
%           .nvar               number of optimization variables
%
%           .gransoFn           function handle to supply to GRANSO.  This 
%                               computes the objective and inequality
%                               constraint, along with their respective
%                               gradients.
%
%           .xToMatrix          Converts vector of variables x to matrix 
%                               M = [ F G; H D ]
%
%           .xToSubMatrices     Converts vector of variables x to matrices
%                               [F,G,H,D] = s.xToSubMatrices(x)
%           
%           .froNormSqFn        Computes the Frobenius norm squared
%                               objective for input matrix M = [ F G; H D ]  
%
%           .specNormSqFn       Computes the spectral norm objective
%                               for input matrix M = [ F G; H D ]  
%   
%           .specAbsRadFn       For input matrix M = [ F G; H D ], this
%                               function computes the spectral
%                               abscissa|radius of matrix F.

    r       = size(X0,1);
    nu      = size(U0,1);
    ny      = size(Y1,1);
  
    Y       = [ X1; Y1 ];
    U       = [ X0; U0 ];
    X1Y1    = [ X1; Y1 ];
    
    % default options
    opts    = struct(   'discrete_time',        true,       ...
                        'frobenius_obj',        true,       ...
                        'regularization',       0,          ...
                        'regularization_type',  1,          ...
                        'require_stable',       true,       ...
                        'stability_margin',     0           );
                    
    % Process user options if any, to determine problem's characteristics                
    if nargin < 6
        user_opts = struct();
    end
    % Define the nested/shared variables that processOpts will set
    [absRadFn,sortFn,normFn,conFn,regFn] = deal([]);
    processOpts();
      
    % shortcut variables, instead of accessing the struct
    discrete            = opts.discrete_time;
    stab_margin         = opts.stability_margin;
    nvar                = (r+nu)*(r+ny);  
    
    % So we can easily replicate all settable parameters in the output
    fields_cell         = fieldnames(opts);
    values_cell         = struct2cell(opts);
    opts_cell           = cell(1,2*length(values_cell));
    opts_cell(1:2:end)  = fields_cell(:);
    opts_cell(2:2:end)  = values_cell(:);
   
    % Output: resuling problem specs and function handles
    prob = struct(                                                      ...
        opts_cell{:},                                                   ...                                                    
        'nvar',                 nvar,                                   ...
        'explicitFn',           @explicit,                              ...
        'gransoFn',             @gransoFn,                              ...
        'xToMatrix',            @toMatrix,                              ...
        'xToSubMatrices',       @(x) toSubMatrices(toMatrix(x)),        ...
        'froNormSqFn',          @frobeniusNormSquared,                  ...
        'specNormFn',           @spectralNorm,                          ...
        'specAbsRadFn',         @spectralAbsRadF                        );
    
    % Nested Functions
                
    function processOpts()
        setUserOption('discrete_time');
        setUserOption('frobenius_obj');
        setUserOption('regularization');
        setUserOption('regularization_type');
        setUserOption('require_stable'); 
        setUserOption('stability_margin');  
        
        if opts.discrete_time
            absRadFn    = @(x) abs(x);
        else
            absRadFn    = @(x) real(x);
        end
        sortFn          = @(x) sort(absRadFn(x),'descend');
        
        if opts.frobenius_obj
            normFn      = @frobeniusNormSquared;
        else
            normFn      = @spectralNorm;
        end
        
        if opts.regularization < 0
            error('opts.regularization must be nonnegative.');
        end
        
        if opts.regularization == 0
            regFn       = @noRegularization;
        else
            reg_amount  = opts.regularization;
            switch opts.regularization_type
                case 1
                    regTermFn = @regSpecNormSquared;
                otherwise
                    error('opts.regularization_type is unrecognized.');
            end
            regFn       = @regularizationFn;
        end
        
        if opts.require_stable
            if opts.stability_margin < 0
                error('opts.stability_margin must be nonnegative.');
            end
            conFn       = @stabilityConstraintFn;
        else
            conFn       = @noConstraint;
        end
        
        function [r,r_grad] = regularizationFn(x)
            [r,r_grad]  = regTermFn(x);
            r           = reg_amount*r;
            r_grad      = reg_amount*r_grad;
        end
    end
            
    function [f,f_grad,ci,ci_grad,ce,ce_grad] = gransoFn(x)
        FGHD            = toMatrix(x);
        
        % Objective value 
        [f,f_grad]      = normFn(FGHD);
        
        % Regularization term for objective, if any
        [reg,reg_grad]  = regFn(x);
        f               = f + reg;
        f_grad          = f_grad + reg_grad;
        
        % Inequality constraint, if any
        [ci,ci_grad]    = conFn(FGHD);
        
        % No equality constraint
        ce              = [];
        ce_grad         = [];
    end

    function [c,c_grad] = stabilityConstraintFn(FGHD)
        [c,c_grad]      = spectralAbsRadF(FGHD);
        c               = c + stab_margin;
        if discrete
            c           = c - 1;
        end
    end
        
    function [F,G,H,D] = explicit()
        FGHD        = (U'\Y');
        [F,G,H,D]   = toSubMatrices(FGHD');
    end

    function [f,g] = frobeniusNormSquared(FGHD)
        Y_MU        = Y - FGHD*U;
        f           = norm(Y_MU,'fro')^2;% + 0.1*norm(FGHD(:),'fro')^2;
        g           = -2*Y_MU*U';
        g           = g(:);
    end

    function [f,g] = spectralNorm(FGHD)
        
        [F,G,H,D]   = toSubMatrices(FGHD);
        FGHD        = [ F*X0 + G*U0
                        H*X0 + D*U0       ];
        [Us,S,Vs]   = svd(FGHD - X1Y1);
        
        % two norm plus associated singular vectors
        f           = S(1);
        u           = Us(:,1);
        v           = Vs(:,1);
        
        % Compute the gradients for each submatrix
        dF          = yhXzGradient(u(1:r),X0*v);
        dG          = yhXzGradient(u(1:r),U0*v);
        dH          = yhXzGradient(u(r+1:end),X0*v);
        dD          = yhXzGradient(u(r+1:end),U0*v);
        
        % Combine them all and then vectorize it to gradient with respect
        % to matrix M
        dM          = [ reshape(dF,r,r)     reshape(dG,r,nu) 
                        reshape(dH,ny,r)    reshape(dD,ny,nu)   ];            
        g           = real(dM(:));
    end
    
    function [f,g] = spectralAbsRadF(FGHD)
        F           = FGHD(1:r,1:r);

        % Compute the right eigenvectors
        [R,dR]      = eig(F,'vector');
        [~,indx]    = sortFn(dR);
        dR          = dR(indx);
        R           = R(:,indx);
        x           = R(:,1);
        z           = dR(1);
        
        % Find the matching left eigenvector
        [L,dL]      = eig(F','vector');
        [~,indx]    = min(abs(dL - conj(z)));
        indx        = indx(1); 
        y           = L(:,indx);
        
        % Spectral abscissa|radius  of F 
        f           = absRadFn(z);       
        g           = eigenvalueGrad(z,x,y);
    end

    function g = eigenvalueGrad(z,x,y)
         % First compute the associated eigenvalue gradient for F only
        [yhx,mu]    = deal(y'*x,1);
        if discrete
            mu      = conj(z)/abs(z);
        end
        gF          = yhXzGradient(y,x)/yhx;
        gF          = real(mu*gF);
        
        % Then place gF in the properly ordered vector for gFGHD
        g           = zeros(r+ny,r+nu);
        g(1:r,1:r)  = reshape(gF,r,r);
        g           = g(:);
    end

    function M = toMatrix(x)     
        M = reshape(x,r+ny,r+nu);
    end

    function [F,G,H,D] = toSubMatrices(M)
        F = M(1:r,1:r);
        G = M(1:r,r+1:end);
        H = M(r+1:end,1:r);
        D = M(r+1:end,r+1:end);  
    end

    function setUserOption(field)
        if isfield(user_opts,field)
            opts.(field) = user_opts.(field);
        end
    end
end

function [c,c_grad] = noConstraint(~)
    c       = [];
    c_grad  = [];
end

function [r,r_grad] = noRegularization(~)
    r       = 0;
    r_grad  = 0;
end

function [r,r_grad] = regSpecNormSquared(x)
    % The input is the following:
    % X = [F G; H D];
    % x = X(:);
    r       = sum(x.^2);
    r_grad  = 2*x;
end