function [halt_log_fn, get_log_fn] = makeHaltLogFunctions(maxit)
%   makeHaltLogFunctions:
%       This is a template class to instruct users how data tracking the
%       progress of granso at every iteration can easily be collected.  
% 
%   USAGE:
%       >> % Set up logging facilities and GRANSO options  
%       >> [halt_log_fn, get_log_fn]    = makeLoggingFunctions(maxit);
%       >> opts.maxit                   = maxit;
%
%       >> % Assign function handle to GRANSO's opts; not a function call!
%       >> opts.halt_log_fn             = halt_log_fn;  
% 
%       >> % Run GRANSO on your problem with opts struct supplied 
%       >> soln = granso(nvar,obj_fn,ineq_fn,eq_fn,opts);
%       
%       >> % Retrieve GRANSO's history of iterates by calling get_log_fn
%       >> granso_log                   = get_log_fn();        
%
%   INPUT:
%       maxit           max number of GRANSO iterations 
%                       Used so one can preallocate arrays here for storing
%                       tracking data at every iterate (to avoid them being
%                       automatically and inefficiently resized while 
%                       while GRANSO is running).
%   
%   OUTPUT:
%       halt_log_fn     function handle to user-defined nested function 
%                       haltLog(), see below, for collecting desired 
%                       tracking data at every iterate of GRANSO
%    
%       get_log_fn      function handle to user-defined nested function
%                       getLog(), see below, for retrieving the collected 
%                       tracking data after GRANSO has terminated.
%
%   See also granso, gransoOptions, and gransoOptionsAdvanced.
% 
%
%   For comments/bug reports, please visit the GRANSO GitLab webpage:
%   https://gitlab.com/timmitchell/GRANSO
%
%   makeHaltLogFunctions.m introduced in GRANSO Version 1.0.
%
% =========================================================================
% |  GRANSO: GRadient-based Algorithm for Non-Smooth Optimization         |
% |  Copyright (C) 2016 Tim Mitchell                                      |
% |                                                                       |
% |  This file is part of GRANSO.                                         |
% |                                                                       |
% |  GRANSO is free software: you can redistribute it and/or modify       |
% |  it under the terms of the GNU Affero General Public License as       |
% |  published by the Free Software Foundation, either version 3 of       |
% |  the License, or (at your option) any later version.                  |
% |                                                                       |
% |  GRANSO is distributed in the hope that it will be useful,            |
% |  but WITHOUT ANY WARRANTY; without even the implied warranty of       |
% |  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        |
% |  GNU Affero General Public License for more details.                  |
% |                                                                       |
% |  You should have received a copy of the GNU Affero General Public     |
% |  License along with this program.  If not, see                        |
% |  <http://www.gnu.org/licenses/>.                                      |
% =========================================================================
    
    % don't change these function handles
    halt_log_fn     = @haltLog;
    get_log_fn      = @getLog;

    index           = 0;
    x_values        = cell(1,maxit);
    f_values        = zeros(1,maxit);
    v_values        = zeros(1,maxit);

    
    % Make your shared variables here to store GRANSO history data
    % EXAMPLE - store history of iterates x_0,x_1,...,x_k

    
    % Only modify the body of logIterate(), not its name or arguments.
    % Store whatever data you wish from the current GRANSO iteration info,
    % given by the input arguments, into shared variables of
    % makeHaltLogFunctions, so that this data can be retrieved after GRANSO
    % has been terminated.
    % 
    % DESCRIPTION OF INPUT ARGUMENTS
    %   iter                current iteration number
    %   x                   current iterate x 
    %   penaltyfn_parts     struct containing the following
    %       OBJECTIVE AND CONSTRAINTS VALUES
    %       .f              objective value at x
    %       .f_grad         objective gradient at x
    %       .ci             inequality constraint at x
    %       .ci_grad        inequality gradient at x
    %       .ce             equality constraint at x
    %       .ce_grad        equality gradient at x
    %       TOTAL VIOLATION VALUES (inf norm, for determining feasibiliy)
    %       .tvi            total violation of inequality constraints at x
    %       .tve            total violation of equality constraints at x
    %       .tv             total violation of all constraints at x
    %       TOTAL VIOLATION VALUES (one norm, for L1 penalty function)
    %       .tvi_l1         total violation of inequality constraints at x
    %       .tvi_l1_grad    its gradient
    %       .tve_l1         total violation of equality constraints at x
    %       .tve_l1_grad    its gradient
    %       .tv_l1          total violation of all constraints at x
    %       .tv_l1_grad     its gradient
    %       PENALTY FUNCTION VALUES 
    %       .p              penalty function value at x
    %       .p_grad         penalty function gradient at x
    %       .mu             current value of the penalty parameter
    %       .feasible_to_tol logical indicating whether x is feasible
    %   d                   search direction
    %   get_BFGS_state_fn   function handle to get the (L)BFGS state data     
    %                       FULL MEMORY: 
    %                       - returns BFGS inverse Hessian approximation 
    %                       LIMITED MEMORY:
    %                       - returns a struct with current L-BFGS state:
    %                           .S          matrix of the BFGS s vectors
    %                           .Y          matrix of the BFGS y vectors
    %                           .rho        row vector of the 1/sty values
    %                           .gamma      H0 scaling factor
    %   H_regularized       regularized version of H 
    %                       [] if no regularization was applied to H
    %   fn_evals            number of function evaluations incurred during
    %                       this iteration
    %   alpha               size of accepted size
    %   n_gradients         number of previous gradients used for computing
    %                       the termination QP
    %   stat_vec            stationarity measure vector                 
    %   stat_val            approximate value of stationarity:
    %                           norm(stat_vec)
    %                       gradients (result of termination QP)
    %   fallback_level      number of strategy needed for a successful step
    %                       to be taken.  See bfgssqpOptionsAdvanced.
    %
    % OUTPUT ARGUMENT
    %   halt                set this to true if you wish optimization to 
    %                       be halted at the current iterate.  This can be 
    %                       used to create a custom termination condition,
    % 
    function halt = haltLog(    iter, x, penaltyfn_parts, d,        ...
                                get_BFGS_state_fn, H_regularized,   ...
                                ls_evals, alpha,                    ...
                                n_gradients, stat_vec, stat_val,    ...
                                fallback_level                      )
          
        % DON'T CHANGE THIS
        % increment the index/count 
        index = index + 1;

        % EXAMPLE:
        % store history of x iterates in a preallocated cell array
        x_values{index} = x;
        f_values(index) = penaltyfn_parts.f;
        v_values(index) = penaltyfn_parts.tv;
        
        % keep this false unless you want to implement a custom termination
        % condition
        halt = false;
        if penaltyfn_parts.tv == 0
            % quit if feasible and it is the first iteration or the
            % objective value is sufficiently close to the original value
            halt = (iter == 0) || (penaltyfn_parts.f < 1000*f_values(1));
        end
        %halt = (penaltyfn_parts.tv == 0) && (iter == 0);
    end

    % Once GRANSO has run, you may call this function to get retreive all
    % the logging data stored in the shared variables, which is populated 
    % by haltLog being called on every iteration of GRANSO.
    function log = getLog()
        % EXAMPLE
        % return x_iterates, trimmed to correct size 
        log.x   = x_values(1:index);
        log.f   = f_values(1:index);
        log.tv  = v_values(1:index);
    end
end