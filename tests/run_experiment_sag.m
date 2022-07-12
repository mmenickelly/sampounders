function run_experiment_sag(macro_seed)

    % EXPERIMENT 10: comparing sag to sam-fo
    
    addpath('~/asynchpounders20/code/sam');
    
    macro_seed = str2num(macro_seed);

    num_micro_seeds = 3;
        
    m = 256; num_epochs = 100;
    
    for type_no = 1:3
    
        if type_no == 1
            experiment_type = 'balanced';
        elseif type_no == 2
            experiment_type = 'progressive';
        elseif type_no == 3
            experiment_type = 'imbalanced';
        end

        for micro_seed = 1:num_micro_seeds
            experiment3(experiment_type,macro_seed,micro_seed,num_epochs,m);
        end
    
    end

end
