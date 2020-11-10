    function alpha = pas_wolfe(f, x, p, alpha0,varargin)
    %The function alpha = pas_wolfe(f, x, p, alpha0,varargin) computes the
    %length of the step of the gradient method that satisfies the
    %Armijo rule (sufficient decrease) and the Wolfe condition
    %(curvature condition).
    %
    %   Input : f      : A function handle to the objective function of the
    %                    problem. f must accept a first input which is of
    %                    the size of x and return as a first output the
    %                    value of the function evaluated at the input and
    %                    as an optional second argument the gradient of the
    %                    function at the input. Additional inputs can be
    %                    accepted by f, those are the ones stored in the
    %                    inputs following 'alpha0'
    %           x      : The current point at which the gradient was
    %                    computed, this point must be a scalar.
    %           p      : The search direction (gradient or something else)
    %                    its size must be the same as x
    %           alpha0 : A starting point for the algorithm.
    %          (option): all additional inputs for 'f' separated by commas,
    %                    this can be used to treat vector function where
    %                    you optimize on one coordinate at a time.   
    %
    %   Output : alpha : The length of the step that satisfies both the
    %                    Armijo and the Wolfe rule(scalar). 
    %
    %Author : François Glineur (?)
    %Imported and commented on the 29/03/2012. 
    %History : 29/03/2012 : The code was modified so that additional
    %                       arguments could be passed to f. 
    %Last Modification : 29/03/2012
    
    % Paramètres (à mettre dans options ?)
    c1 = 1e-4;
    c2 = 0.9;
    theta = 1.5;

    alpha = alpha0;             % alpha_i (alpha actuel)
    alpha_prec = 0;             % alpha_{i-1} (alpha précédent)
    [phi_val0, phi_prime0] = phi(0);
    if phi_prime0 >= 0
        warning('WOLFE:WEIRD','Avertissement: pas de Wolfe cherché avec un phi''(0) qui n''est pas négatif');
        alpha = alpha0;
        return;
    end
    phi_valprec = inf;          % astuce pour ne jamais déclencher à la première itération la seconde condition du or 
    while 1
        [phi_val, phi_prime] = phi(alpha);
        if phi_val > phi_val0 + alpha*c1*phi_prime0 || phi_val > phi_valprec
            alpha = zoom_wolfe(alpha_prec, alpha);
            return;
        elseif abs(phi_prime) <= -c2*phi_prime0
            return;
        elseif phi_prime >= 0
            alpha = zoom_wolfe(alpha, alpha_prec);
            return,
        else
            alpha_prec = alpha;
            phi_valprec = phi_val;
            alpha = alpha*theta;
        end
    end

        function alpha = zoom_wolfe(alpha_l, alpha_h)
            while 1
                alpha = (alpha_l+alpha_h)/2;        % Améliorer par interpolation ?
                [phi_val, phi_prime] = phi(alpha);
                if phi_val > phi_val0 + alpha*c1*phi_prime0 || phi_val > phi(alpha_l)
                    alpha_h = alpha;
                elseif abs(phi_prime) <= -c2*phi_prime0 || abs(phi_prime) < 1e-7 % !! nécessaire sinon parfois qd phi_prime0 est trop petit impossible à obtenir (erreurs numériques ?)
                    alpha
                    return;
                elseif phi_prime*(alpha_h-alpha_l) >= 0
                    alpha_h = alpha_l;
                    alpha_l = alpha;
                else
                    alpha_l = alpha;
                end
            end
        end

        function [phi_val, phi_prime] = phi(alpha)
            if nargout == 1
                phi_val = f(x + alpha*p,varargin{:});
            else
                [phi_val, grad] = f(x + alpha*p,varargin{:});
                phi_prime = grad'*p;
            end
        end
    end

