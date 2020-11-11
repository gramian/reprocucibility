function RUNME(id)

    myWorks = {'HimO13', ...      % A Unified Software Framework for Empirical Gramians
               'HimO14', ...      % Cross-Gramian-Based Combined State and Parameter Reduction for Large-Scale Control Systems
               'HimO14a', ...     % Model Reduction for Complex Hyperbolic Networks
               'HimO14b', ...     % Combined State and Parameter Reduction
               'HimO15', ...      % The Empirical Cross Gramian for Parametrized Nonlinear Systems
               ...
               'HimO15b', ...     % Data-driven combined state and parameter reduction for inverse problems
               'HimO16', ...      % A note on the cross Gramian for non-symmetric systems
               'FehHHetal16', ... % Best Practices for Replicability, Reproducibility and Reusability of Computer-Based Experiments Exemplified by Model Reduction Software
               ...
               'HimO17', ...      % Cross-Gramian-Based Model Reduction: A Comparison
               ...
               'HimLRetal17', ... % Fast Low-Rank Empirical Cross Gramians
               ...
               'Him18', ...       % emgr - Empirical Gramian Framework
               ...
               'BenHM18', ...     % On Reduced Input-Output Dynamic Mode Decomposition
               ...
               'BenH19',
               ... % tbc
              };

    if isequal(id,'all')

       for m = myWorks

           RUNME(m{:});
       end%for

       return;

    elseif not(any(strcmp(id,myWorks)));

        fprintf('Unknown id, use id=\n\n');
        for m = myWorks

            fprintf(' %s\n',m{1});
        end%for
        fprintf(' all\n\n');
        return
    end%if

    addpath(id);

    switch(id)

        case 'HimO13'
            unified(2);

        case 'HimO14'
            mpe();
            mpe_plots;

        case 'HimO14a'
            hypnet();

        case 'HimO14b' 
            gamm14mri(5);
            gamm14eeg(3);

        case 'HimO15'
            mathmod();
            plot_mathmod;

        case 'HimO15b'
            acom();
            out_acom;

        case 'HimO16'
            nonsym(0);
            nonsym(1);
            nonsym(2);
            nonsym(3);
            nonsym(4);
            nonsym(5);

        case 'FehHHetal16'
            cd('FehHHetal16');
            runme = @RUNME;
            cd('..');
            runme();

        case 'HimO17'
            cd('HimO17');
            runme = @RUNME;
            cd('..');
            runme();

        case 'HimLRetal17'
            cd('HimLRetal17');
            runme = @RUNME;
            plotme = @PLOTME;
            cd('..');
            runme();
            plotme();

        case 'Him18'
            cd('Him18');
            runme = @RUNME;
            cd('..');
            runme();

        case 'BenHM18'
            cd('BenHM18');
            addpath('GRANSO')
            runme = @RUNME;
            cd('..');
            runme(0);
            runme(1);

        case 'BenH19'
            cd('BenH19');
            runme = @RUNME;
            cd('..');
            runme();

    end%switch

    rmpath(id);
end
