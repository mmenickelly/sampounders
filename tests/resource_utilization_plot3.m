function [bhl,medtable] = resource_utilization_plot3(which_solvers,num_macro_seeds,num_micro_seeds,num_epochs,m,test_function,experiment_type,solver,convergence,batchsizes,resource_sizes,gate)

%% INPUTS:
% H - nf x np x ns from run

num_batchsizes = length(batchsizes);
num_probs = num_macro_seeds*num_micro_seeds;

% only doing this now to get prob_min: 
if strcmp(solver,'fo')
    masterHf = zeros(num_probs,2*num_batchsizes-1,num_epochs*m);
    masterHg = zeros(num_probs,2*num_batchsizes-1,num_epochs*m);
else
    masterHf = zeros(num_probs,2*num_batchsizes-1,num_epochs*m*(m+1));
    masterHg = zeros(num_probs,2*num_batchsizes-1,num_epochs*m*(m+1));
end

for macro_seed = 1:num_macro_seeds
    for micro_seed = 1:num_micro_seeds
        filename = strcat('results/',test_function,'_',experiment_type,'_',solver,'_',num2str(macro_seed),'_',num2str(micro_seed),'_',num2str(m),'.mat');
        load(filename);
        masterHf(num_micro_seeds*(macro_seed-1) + micro_seed,:,:) = Hf;
        masterHg(num_micro_seeds*(macro_seed-1) + micro_seed,:,:) = Hg;
    end
end

if strcmp(convergence,'fval')
    masterHf = permute(masterHf,[3 1 2]);
    prob_min = min(min(masterHf),[],3);
elseif strcmp(convergence,'gval')
    masterHg = permute(masterHg,[3 1 2]);
    prob_min = min(min(masterHg),[],3);
end  
% end prob_min acquisition

num_resource_sizes = length(resource_sizes);
[nf,np] = size(Evals{1}); 

%Plotting commands:
colors=get(gca,'colororder');
lines={'-','-.','--'};
markers={'s','o','^','v','p','<','x','h','+','d','*'};
LW=2;
FS=16;
MS=8;

target_values = prob_min + gate;

resourcelabels = {};
for j = 1:num_resource_sizes
    resourcelabels{j} = num2str(batchsizes(j));
end

% medtable
medtable = zeros(num_resource_sizes);
resource_use = nan(np,num_resource_sizes); % matrix to be populated
for s = 1:num_resource_sizes 
    p = 0;
    for macro_seed = 1:num_macro_seeds
        for micro_seed = 1:num_micro_seeds
            
            p = p + 1;
            
            filename = strcat('results/',test_function,'_',experiment_type,'_',solver,'_',num2str(macro_seed),'_',num2str(micro_seed),'_',num2str(m),'.mat');
            load(filename);
    
            if strcmp(convergence,'fval')
                Hf = permute(Hf,[2 1]);
                H = Hf;
            elseif strcmp(convergence,'gval')
                Hg = permute(Hg,[2 1]);
                H = Hg; 
            end  
        
            data = H(:,which_solvers(s));
            Eval = Evals{which_solvers(s)};
            
            under_target = find(data < target_values(p),1);
            if isempty(under_target)
                under_target = nf;
            end
            % SHOULD COUNT COMPUTATIONAL ROUNDS NEEDED TO REACH under_target
            a = sum(Eval,2); 
            cumsuma = cumsum(a);
            ind = find(cumsuma > under_target,1);
            if isempty(ind)
                ind = nf;
            end
            rounds = sum(ceil(a(1:ind)/resource_sizes(s)))*ceil((resource_sizes(s)./batchsizes));
            resource_use(p,:) = rounds;
        
        
        end
    end

    % CREATE A CSV TO MAKE INTO A TABLE
    med = quantile(resource_use,0.5);
    medtable(s,:) = med;
    
    resource_use = log2(resource_use); % why log2 and not some other log base? have to remember 

    med = quantile(resource_use,0.5);
    bhl{s} = semilogx(batchsizes,med,'Color',colors(mod(s-1,7)+1,:),'Marker',markers{mod(s-1,11)+1},'MarkerSize',MS,'LineStyle',lines{mod(s-1,3)+1},'LineWidth',LW);
    ax = gca;
    ax.XTick = batchsizes;
    ax.XTickLabels = resourcelabels; 
    hold on
end

set(gca,'FontSize',FS)
set(gca,'XLim',[batchsizes(1),batchsizes(end)]); 
% xlims = get(gca,'XLim'); ylims = get(gca,'YLim');
% xlims(1) = xlims(1) - .25; xlims(2) = xlims(2) + .25;
% ylims(1) = ylims(1) - 1; ylims(2) = ylims(2) + 1;
% set(gca,'XLim',xlims); set(gca,'YLim',ylims);
% axis tight

end
