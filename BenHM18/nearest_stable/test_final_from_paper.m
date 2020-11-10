% Grand test final : this file can be run to make exactly* the same
% experiment as in our paper 
% 
% *except for the random type matrices that are regenerated at each run 

clear all; clc; 

addpath(genpath('.\Orbandexivry - SuccConv'));
addpath(genpath('.\Overton - BFGS')); 

maxiter = 1e16; 

dimn   = [10 20 50 100];   % dimensions
timmax = [20 100 300 600]; % time limits 
% To have a quick test, use 
% timmax = [2 2 2 2];

numtm = 4;             % number of types of matrices
numalg = 5;            % number of tested algorithms 

disp('Estimated total time in hours') 
length(dimn) * numtm * numalg * mean(timmax)  / 3600 

% loop over different size of matrices 
for dn = 1 : length(dimn) 
    n = dimn(dn); 
    % loop over the different types of matrices 
    for tm = 1 : numtm
        if tm == 1
            A = fxorban(n,-0.1);
        elseif tm == 2
            A = grcar(n,3); 
        elseif tm == 3
            A = randn(n); 
        elseif tm == 4
            A = rand(n); 
        end

        All{dn}{tm} = A; 
        J = (A-A')/2; 
        R = projectPSD(J-A); 
        e0(dn,tm) = norm(A-(J-R),'fro');  
        
        % loop over the different alorithms
        for alg = 1 : numalg
            clear X t e 
            if alg == 4 && n > 50
                X = []; 
                e = e0(dn,tm); 
                t = 0; 
            else
                [X,e,t] = neareststablealgo( A, maxiter, timmax(dn), alg  ); 
            end
            Xall{dn}{tm}{alg} = X; 
            eall{dn}{tm}{alg} = e; 
            tall{dn}{tm}{alg} = t; 
            eend( dn, tm, alg  ) = min(e);  
        end
    end
end 
save results_rand


% Tables 
for dn = 1 : length(dimn) 
    n = dimn(dn); 
    for tm = 1 : numtm
                % Display results for each matrix and each dimension
        fprintf('----------- type =%2.0f, n = %3.0f  ----------- \n', tm, n)
        fprintf('  Algorithm  |   Error  |  Iterations \n')
        fprintf('---------------------------------------------- \n')
        fprintf(' Init. Error |   %2.2f   | %1.0f \n',  e0(dn,tm) , 0 );
        fprintf(' BCD         |   %2.2f   | %1.0f \n',  eend(dn, tm, 1), length(eall{dn}{tm}{1}) );
        fprintf(' Grad       |   %2.2f   | %1.0f \n',  eend(dn, tm, 2), length(eall{dn}{tm}{2}) );
        fprintf(' FGM         |   %2.2f   | %1.0f \n',  eend(dn, tm, 3), length(eall{dn}{tm}{3}) );
        fprintf(' SuccConv    |   %2.2f   | %1.0f \n',  eend(dn, tm, 4), length(eall{dn}{tm}{4}) );
        fprintf(' BFGS        |   %2.2f   | %1.0f \n',  eend(dn, tm, 5), length(eall{dn}{tm}{5}) ); 
        fprintf('---------------------------------------------- \n')
    end
end

% Figures 
set(0, 'DefaultAxesFontSize', 18);
set(0, 'DefaultLineLineWidth', 2);
for dn = 1 : length(dimn) 
    n = dimn(dn); 
    
    figure; 
    for tm = 1 : numtm
        subplot( 2, 2 , tm ); 
        
        % FiX
        set(0, 'DefaultLineLineWidth', 1);
        plot( tall{dn}{tm}{4}, eall{dn}{tm}{4} ,'ko--', 'MarkerSize', 5); hold on;   
        % Overton
        set(0, 'DefaultLineLineWidth', 5);
        plot( tall{dn}{tm}{5}, eall{dn}{tm}{5} ,'-', 'color',[0.4,0.4,0.4]); 
        % BCD
        set(0, 'DefaultLineLineWidth', 1);
    	plot( tall{dn}{tm}{1}, eall{dn}{tm}{1} ,'bx--', 'MarkerSize', 8); 
        % Grad
        set(0, 'DefaultLineLineWidth', 2);
        plot( tall{dn}{tm}{2}, eall{dn}{tm}{2} ,'r'); 
         % FG 
        plot( tall{dn}{tm}{3}, eall{dn}{tm}{3} ,'r--');
        
        xlabel('Time (s.)'); 
        ylabel('Error'); 
        axis([0 timmax(dn) min(eend(dn,tm,:))*0.975 e0(dn,tm)*1.025])
        % title
        tit = sprintf('type = %1.0f, n = %1.0f',tm,n); 
        title(tit); 
        if tm == 1
            legend('SuccConv','BFGS','BCD','Grad','FGM'); 
        end
    end
    %set(gca,'position',[0 0 1 1],'units','normalized') 
    iptsetpref('ImshowBorder','tight');
   
end

% % LaTeX Tables 
% for dn = 1 : length(dimn) 
%     n = dimn(dn); 
%     %for tm = 1 : numtm
%         fprintf('\\begin{center}  \n \\begin{table}[h!] \n \\begin{center} \n'); 
%         fprintf('\\caption{Comparison of the algorithms for matrices $$A$$ of size %1.0f. The table displays the final error obtained by each algorithm and, in brackets, the number of iterations performed.} \n', n); 
%         fprintf('\\label{results%1.0f%1.0f} \n \\begin{tabular}{|c|c|c|c|c|} \n \\hline', tm, n); 
%         fprintf('    &  Type 1 & Type 2 & Type 3 & Type 4    \\\\ \n \\hline \n')
%         fprintf(' Initial error &   %2.2f & %2.2f  & %2.2f  & %2.2f   \\\\ \n',  e0(dn,1) ,  e0(dn,2),  e0(dn,3),  e0(dn,4) );
%         fprintf(' BCD         &   %2.2f (%1.0f) &   %2.2f (%1.0f) &   %2.2f (%1.0f) &   %2.2f (%1.0f)  \\\\ \n', eend(dn, 1, 1), length(eall{dn}{1}{1}), eend(dn, 2, 1), length(eall{dn}{2}{1}), eend(dn, 3, 1), length(eall{dn}{3}{1}), eend(dn, 4, 1), length(eall{dn}{4}{1}) );
%         fprintf(' Grad        &   %2.2f (%1.0f) &   %2.2f (%1.0f) &   %2.2f (%1.0f) &   %2.2f (%1.0f) \\\\ \n', eend(dn, 1, 2), length(eall{dn}{1}{2}), eend(dn, 2, 2), length(eall{dn}{2}{2}), eend(dn, 3, 2), length(eall{dn}{3}{2}), eend(dn, 4, 2), length(eall{dn}{4}{2}) );
%         fprintf(' Fast Grad   &   %2.2f (%1.0f) &   %2.2f (%1.0f) &   %2.2f (%1.0f) &   %2.2f (%1.0f) \\\\ \n', eend(dn, 1, 3), length(eall{dn}{1}{3}), eend(dn, 2, 3), length(eall{dn}{2}{3}), eend(dn, 3, 3), length(eall{dn}{3}{3}), eend(dn, 4, 3), length(eall{dn}{4}{3}) );
%         fprintf(' SuccConv    &   %2.2f (%1.0f) &   %2.2f (%1.0f) &   %2.2f (%1.0f) &   %2.2f (%1.0f) \\\\ \n', eend(dn, 1, 4), length(eall{dn}{1}{4}), eend(dn, 2, 4), length(eall{dn}{2}{4}), eend(dn, 3, 4), length(eall{dn}{3}{4}), eend(dn, 4, 4), length(eall{dn}{4}{4}) );
%         fprintf(' BFGS        &   %2.2f (%1.0f) &   %2.2f (%1.0f) &   %2.2f (%1.0f) &   %2.2f (%1.0f) \\\\ \\hline \n', eend(dn, 1, 5), length(eall{dn}{1}{5}), eend(dn, 2, 5), length(eall{dn}{2}{5}), eend(dn, 3, 5), length(eall{dn}{3}{5}), eend(dn, 4, 5), length(eall{dn}{4}{5}) );
%         fprintf('\\end{tabular} \n \\end{center} \n \\end{table} \n \\end{center} \n'); 
%     %end
% end