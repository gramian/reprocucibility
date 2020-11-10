% Algorithms for finding the nearest stable matrix

function [X,e,t] = neareststablealgo(A,maxiter,timmax,algo)

if algo == 1 % BCD
    [J,R,Q,e,t] = stableLinearBCD(A,maxiter,timmax); 
    X = (J-R)*Q; 
elseif algo == 2 % Grad
    [J,R,Q,e,t] = stableLinearGrad(A,maxiter,timmax); 
    X = (J-R)*Q; 
elseif algo == 3 % FGM 
    [J,R,Q,e,t] = stableLinearFGM(A,maxiter,timmax); 
    X = (J-R)*Q; 
elseif algo == 4 % SuccConv
    [X,e,t]= StableMatMain(A,maxiter,timmax); 
elseif algo == 5 % BFGS 
    [X,e,t] = nearstabdemo(A, maxiter, timmax); 
end