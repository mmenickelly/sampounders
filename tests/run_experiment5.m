function run_experiment5(macro_seed)

    % EXPERIMENT 2: POUNDERS on rosenbrock
    
    addpath('~/sampounders/');

    macro_seed = str2num(macro_seed);

    num_micro_seeds = 3;
      
    batchsizes = [1,2,4,8,16];
        
    m = max(batchsizes); 
    
    num_epochs = 100;
   
    test_function = 'rosenbrock';

    solver = 'pounders';
    
    for type_no = 1:3
    
        if type_no == 1
            experiment_type = 'balanced';
        elseif type_no == 2
            experiment_type = 'progressive';
        elseif type_no == 3
            experiment_type = 'imbalanced';
        end

        for micro_seed = 1:num_micro_seeds
            experiment2(test_function,experiment_type,solver,macro_seed,micro_seed,batchsizes,num_epochs,m);
        end
    
    end
end
