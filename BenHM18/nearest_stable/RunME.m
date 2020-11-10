% Comparing algorithms for nearest stable matrix 
clear all; clc; 

addpath(genpath('.\Orbandexivry - SuccConv'));
addpath(genpath('.\Overton - BFGS'));

% Maximum number of iteration for each algorithm
maxiter = 1e16; 
% Maximum cpu time for each algorithm
timmax = 2; 
% dimension of the input matrix
n = 10; 
% Input matrix: you can choose on which type of matrices you want the code
% to run on (cf. paper for more details). 
typemat = 2; 
if typemat == 1 
    A = fxorban(n,-0.1);
elseif typemat == 2
    A = grcar(n,3); 
elseif typemat == 3
    A = randn(n); 
elseif typemat == 4
    A = rand(n); 
end

% Initialization 
J = (A-A')/2; 
R = projectPSD(J-A); 
e0 = norm(A-(J-R),'fro');  
        
% Run the 5 alorithms 
algor{1} = 'BCD'; 
algor{2} = 'Grad'; 
algor{3} = 'FGM';
algor{4} = 'SuccConv';
algor{5} = 'BFGS'; 
for alg = 1 : 5
    if alg == 4 && n > 50
        X = [];
        e = e0(dn,tm);
        t = 0;
    else
        disp('******************************************'); 
        fprintf(['Running ' algor{alg} ' ... \n']); 
        [X,e,t] = neareststablealgo( A, maxiter, timmax , alg  );
        fprintf(['Done with ' algor{alg} '.\n']); 
        disp('******************************************'); 
    end
    Xall{alg} = X;
    eall{alg} = e;
    tall{alg} = t;
    eend( alg  ) = min(e);
end

% Table with the smallest objective function for each algoirthm 
disp('******************************************'); 
fprintf('Table with final errors and number of iterations, for n = %3.0f: \n', n)
fprintf('--------------------------------------- \n')
fprintf('  Algorithm  |   Error  |  Iterations \n')
fprintf('--------------------------------------- \n')
fprintf(' Init. Error |   %2.2f   | %1.0f \n',  e0 , 0 );
fprintf(' BCD         |   %2.2f   | %1.0f \n',  eend(1), length(eall{1}) );
fprintf(' Grad        |   %2.2f   | %1.0f \n',  eend(2), length(eall{2}) );
fprintf(' FGM         |   %2.2f   | %1.0f \n',  eend(3), length(eall{3}) );
fprintf(' SuccConv    |   %2.2f   | %1.0f \n',  eend(4), length(eall{4}) );
fprintf(' BFGS        |   %2.2f   | %1.0f \n',  eend(5), length(eall{5}) );
fprintf('--------------------------------------- \n')

% Figure displaying the evolution of the objective function 
set(0, 'DefaultAxesFontSize', 18);
set(0, 'DefaultLineLineWidth', 2);
figure;
% FiX Orban
set(0, 'DefaultLineLineWidth', 1);
plot( tall{4}, eall{4} ,'ko--', 'MarkerSize', 5); hold on;
% Overton
set(0, 'DefaultLineLineWidth', 5);
plot( tall{5}, eall{5} ,'-', 'color',[0.4,0.4,0.4]);
% BCD
set(0, 'DefaultLineLineWidth', 1);
plot( tall{1}, eall{1} ,'bx--', 'MarkerSize', 8);
% Grad
set(0, 'DefaultLineLineWidth', 2);
plot( tall{2}, eall{2} ,'r');
% FG
plot( tall{3}, eall{3} ,'r--');
xlabel('Time (s.)');
ylabel('Error');
axis([0 timmax min(eend(:))*0.975 e0*1.025])
% title
tit = sprintf('n = %1.0f',n);
title(tit);
legend('SuccConv','BFGS','BCD','Grad','FGM'); 