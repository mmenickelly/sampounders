function run_experiment4(macro_seed)

    % EXPERIMENT 1: FO on logistic
    
    addpath('~/sampounders/');
    
    macro_seed = str2num(macro_seed);

    num_micro_seeds = 3;
      
    batchsizes = [1,2,4,8,16,32,64,128,256];
        
    m = 256; num_epochs = 100;
   
    test_function = 'logistic';

    solver = 'fo';
    
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
