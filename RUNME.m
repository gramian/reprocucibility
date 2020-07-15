function RUNME(id)

    myWorks = {'HimO13', ...      % A Unified Software Framework for Empirical Gramians
               'HimO14', ...      % Cross-Gramian-Based Combined State and Parameter Reduction for Large-Scale Control Systems
               'HimO14a', ...     % Model Reduction for Complex Hyperbolic Networks
               'HimO14b', ...     % Combined State and Parameter Reduction
               'HimO15', ...      % The Empirical Cross Gramian for Parametrized Nonlinear Systems
               ... % tbc
              };

    if isequal(id,'all')

       for m = myWorks

           RUNME(m);
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

    end%switch

    rmpath(id);
end
