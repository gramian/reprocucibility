function [X,diffs] = getNearestStable(A,user_opts)

    opts.maxit      = 1000;
    opts.max_time   = inf;
    opts.method     = 3;        % supposed to be best
    if nargin > 1 
        setUserOption('maxit');
        setUserOption('max_time');
        setUserOption('method');
    end
    if opts.maxit <= 0 || mod(opts.maxit,1) > 0
        error('opts.maxit must be a positive integer.');
    end
    if opts.max_time <= 0
        error('opts.max_time must be a positive value.');
    end
    if mod(opts.method,1) > 0 || opts.method < 1 || opts.method > 3
        error('opts.method must be in {1,2,3}.');
    end
    
    switch opts.method
        case 1 
            getStableFn = @stableLinearBCD;
        case 2
            getStableFn = @stableLinearGrad;
        case 3
            getStableFn = @stableLinearFGM;
    end
    
    [J,R,Q,diffs]   = getStableFn(A,opts.maxit,opts.max_time);
    X               = (J-R)*Q; 
    
    function setUserOption(field)
        if isfield(user_opts,field)
            opts.(field) = user_opts.(field);
        end
    end
end