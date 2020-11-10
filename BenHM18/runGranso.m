function [soln,t] = runGranso(problem,granso_opts)
%   Takes the the output prob = makeProblem(...) and runs GRANSO on this
%   problem and returns GRANSO's soln output argument.
% 
%   GRANSO options can be provided by the optional second argument.

    if nargin < 2 || ~isstruct(granso_opts) 
        granso_opts = [];
    end
    granso_opts = gransoOptions(problem.nvar,granso_opts);
    granso_opts.quadprog_info_msg = false;
    
    t = tic();
    soln = granso(problem.nvar,problem.gransoFn,granso_opts);
    t = toc(t);
 
end