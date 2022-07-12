function run_experiment3(macro_seed)

    % EXPERIMENT 3: POUNDERS on cube
    
    addpath('~/asynchpounders20/code/sam');
    addpath('~/asynchpounders20/code/sam/minq5');

    macro_seed = str2num(macro_seed);

    num_micro_seeds = 3;
      
    batchsizes = [1,2,4,8,16];
        
    m = max(batchsizes); 
    
    num_epochs = 100;
   
    test_function = 'cube';

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
            experiment4(test_function,experiment_type,solver,macro_seed,micro_seed,batchsizes,num_epochs,m);
        end
    
    end
    % if plotting a data profile as part of the run of this script
%     if print_fig
%         tau = 1e-5; logplot = 1;
%         num_batchsizes = length(batchsizes);
%         %num_nus = length(nus); 
%         num_probs = num_macro_seeds*num_micro_seeds;
% 
%         masterHf = zeros(num_probs,2*num_batchsizes-1,kappa*m);
%         masterHg = zeros(num_probs,2*num_batchsizes-1,kappa*m);
%         for macro_seed = 1:num_macro_seeds
%             for micro_seed = 1:num_micro_seeds
%                 filename = strcat(test_function,'_',experiment_type,'_',num2str(macro_seed),'_',num2str(micro_seed),'.mat');
%                 load(filename);
%                 masterHf(num_micro_seeds*(macro_seed-1) + micro_seed,:,:) = Hf(1:2*num_batchsizes-1,:);
%                 masterHg(num_micro_seeds*(macro_seed-1) + micro_seed,:,:) = Hg(1:2*num_batchsizes-1,:);
%             end
%         end
%         masterHf = permute(masterHf,[3 1 2]);
%         masterHg = permute(masterHg,[3 1 2]);
% 
%         prob_min = min(min(masterHf),[],3); 
% 
%         figure(1)
%         if strcmp(experiment_type,'imbalanced') || strcmp(experiment_type,'progressive') || strcmp(experiment_type,'balanced') 
%             % component functions
%             budget_units = m*ones(num_probs,1);
%         else
%             % simplex gradients x component functions
%             budget_units = m*(m+1)*ones(num_probs,1);
%         end
%         [hlf_uni,hlfs_uni] = data_profile(masterHf(:,:,1:num_batchsizes),budget_units,prob_min,tau,logplot);
%         legend(hlfs_uni,label_strings(1:num_batchsizes),'Location','NorthWest');
%         xl = xlim;
% 
%         figure(2)
%         [hlf_arb,hlfs_arb] = data_profile(masterHf(:,:,[num_batchsizes+1:2*num_batchsizes-1,num_batchsizes]),budget_units,prob_min,tau,logplot);
%         legend(hlfs_arb,label_strings([num_batchsizes+1:2*num_batchsizes-1,num_batchsizes]),'Location','NorthWest');
%         xl2 = xlim;
%         xlmin = min(xl(1),xl2(1));
%         xlmax = max(xl(1),xl2(2));
%         xlim([xlmin xlmax])
%         savestr = strcat(test_function,'_',experiment_type,'_','arbitrary_f.eps');
%         saveas(2,savestr,'epsc')
%         
%         figure(1)
%         xlim([xlmin xlmax])
%         if strcmp(experiment_type,'imbalanced') || strcmp(experiment_type,'progressive') || strcmp(experiment_type,'balanced') 
%             % component functions
%             xlabel('(component function/gradient evaluations)/$p$','interpreter','latex')
%         else
%             % simplex gradients x component functions
%             xlabel('(component function evaluations)/$(p\times(n+1))$','interpreter','latex');
%         end
%         savestr = strcat(test_function,'_',experiment_type,'_','uniform_f.eps');
%         saveas(1,savestr,'epsc')
%         
%         prob_min = min(min(masterHg),[],3); 
% 
%         figure(3)
%         [hlg_uni,hlgs_uni] = data_profile(masterHg(:,:,1:num_batchsizes),budget_units,prob_min,tau,logplot);
%         legend(hlgs_uni,label_strings(1:num_batchsizes),'Location','NorthWest');
%         xl = xlim;
%         
%         figure(4)
%         [hlg_arb,hlgs_arb] = data_profile(masterHg(:,:,[num_batchsizes+1:2*num_batchsizes-1,num_batchsizes]),budget_units,prob_min,tau,logplot);
%         legend(hlgs_arb,label_strings([num_batchsizes+1:2*num_batchsizes-1,num_batchsizes]),'Location','NorthWest');
%         xl2 = xlim;
%         xlmin = min(xl(1),xl2(1));
%         xlmax = max(xl(1),xl2(2));
%         xlim([xlmin xlmax])
%         savestr = strcat(test_function,'_',experiment_type,'_','arbitrary_g.eps');
%         saveas(4,savestr,'epsc')
%         
%         figure(3)
%         xlim([xlmin xlmax])
%         if strcmp(experiment_type,'imbalanced') || strcmp(experiment_type,'progressive') || strcmp(experiment_type,'balanced') 
%             % component functions
%             xlabel('(component function/gradient evaluations)/$p$','interpreter','latex')
%         else
%             % simplex gradients x component functions
%             xlabel('(component function evaluations)/$(p\times(n+1))$','interpreter','latex');
%         end
%         savestr = strcat(test_function,'_',experiment_type,'_','uniform_g.eps');
%         saveas(3,savestr,'epsc')
%     end
end
