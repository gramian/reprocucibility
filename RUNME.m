function RUNME(id)

    myWorks = {'HimO13', ...      % A Unified Software Framework for Empirical Gramians
               'HimO14', ...      % Cross-Gramian-Based Combined State and Parameter Reduction for Large-Scale Control Systems
               'HimO14a', ...     % Model Reduction for Complex Hyperbolic Networks
               ... % tbc
              };

    if isequal(id,'all')

       for m = myWorks

           RUNME(m);
       end%for

       return;

    elseif all(not(strcmp(id,myWorks)));

        fprintf('Unknown id, use id=\n\n');
        for m = myWorks

            fprintf(' %s\n',m);
        end%for
        fprintf(' all\n\n');
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

    end%switch

    rmpath(id);
end
