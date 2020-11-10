function g = yhXzGradient(y,z)
% get the gradient of y'*X*z with respect to a matrix X and vectors y and z
    g = conj(y)*z.';
    g = g(:);
end