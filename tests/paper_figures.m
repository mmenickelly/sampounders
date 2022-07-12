function paper_figures(which_plots)

if nargin == 0
    which_plots = 1:6;
end

% THESE ARE THE CORRECT SETTINGS FOR ALL PLOTS:
types = {'balanced','progressive','imbalanced'};
num_macro_seeds = 30;
num_micro_seeds = 3;
num_epochs = 100;
convergence = 'fval';
num_probs = num_macro_seeds*num_micro_seeds; 
rup_title = '$\tau =  $';

gates = logspace(-6,-1,100);
to_label = [-6,-5,-4,-3,-2,-1];

rup_gates = [1e-1,1e-3,1e-5,1e-7];
FS = 16;

mode = 'eval';

%% FIGURE 1 PLOTS
if ismember(1,which_plots)
% SPECIFIC TO EXPERIMENT 1: 
m = 256;
%batchsizes = [1,2,4,8,16,32,64,128,256];
batchsizes = [1,4,16,64];
uniform_solvers = [1,3,5,7]; dynamic_solvers = 9 + uniform_solvers;
ns = length(uniform_solvers); 
test_function = 'logistic';
solver = 'fo';
titles = {'Balanced datasets', 'Progressive datasets', 'Imbalanced datasets'};

legend_str = {};
for j = 1:ns
   legend_str{2*j-1} = strcat('Uniform, $r=$',num2str(batchsizes(j)));
   legend_str{2*j} = strcat('Dynamic, $r=$',num2str(batchsizes(j)));
end

width = 25;

figure; hlt = tiledlayout(1,3);
for j=1:3
    experiment_type = types{j};
    masterHf = zeros(num_probs,2*ns,num_epochs*m);
    masterHg = zeros(num_probs,2*ns,num_epochs*m);
    for macro_seed = 1:num_macro_seeds
        for micro_seed = 1:num_micro_seeds
            filename = strcat('results/',test_function,'_',experiment_type,'_',solver,'_',num2str(macro_seed),'_',num2str(micro_seed),'_',num2str(m),'.mat');
            load(filename);
            masterHf(num_micro_seeds*(macro_seed-1) + micro_seed,1:ns,:) = Hf(uniform_solvers,:);
            masterHf(num_micro_seeds*(macro_seed-1) + micro_seed,ns+1:2*ns,:) = Hf(dynamic_solvers,:);
            masterHg(num_micro_seeds*(macro_seed-1) + micro_seed,1:ns,:) = Hg(uniform_solvers,:);
            masterHg(num_micro_seeds*(macro_seed-1) + micro_seed,ns+1:2*ns,:) = Hg(dynamic_solvers,:);
        end
    end
    
    if strcmp(convergence,'fval')
        masterHf = permute(masterHf,[3 1 2]);
        prob_min = min(min(masterHf),[],3);
        H = masterHf;
    elseif strcmp(convergence,'gval')
        masterHg = permute(masterHg,[3 1 2]);
        prob_min = min(min(masterHg),[],3);
        H = masterHg; 
    end
    nexttile;
    [hlmu,hlmd] = plot_trajectory_percentiles2(H,prob_min,1:ns,ns+1:2*ns,width,mode,batchsizes,m);
    hold on;
    if j == 1
        hlm = [];
        for s = 1:ns
            hlm = cat(1,hlm,hlmu(s)); hlm = cat(1,hlm,hlmd(s));
        end
        legend(hlm,legend_str,'interpreter','latex','FontSize',FS);
    end
    if j == 1
        ylabel('$f(\mathbf{x}^k) - f(\mathbf{x}^*)$','interpreter','latex','FontSize',FS);
    end
    if j == 2
        if strcmp(mode,'eval')
            xlabel('Effective data passes','FontSize',FS);
        elseif strcmp(mode,'batch')
            xlabel('Batches of $r$ evaluations','FontSize',FS,'interpreter','latex');
        end
    end
    if j == 1
        ylim([1e-17,1e5])
    elseif j == 2
        ylim([1e-10,1e2])
    elseif j == 3
        ylim([1e-16,1e2])
    end
    xlim([0,100])
    title(titles{j},'FontSize',FS);
 
end
hlt.Padding = 'compact'; hlt.TileSpacing = 'compact'; 
pause() % MAKE IT CAMERA READY DURING THE PAUSE
saveimgname = strcat('images/logistic_uni_v_dyn_',test_function,'_','.eps');
saveas(gcf,saveimgname,'epsc');
end

%% EXPERIMENT 2 PLOTS
if ismember(2,which_plots)
% SPECIFIC TO EXPERIMENT 2: 
m = 16;
%batchsizes = [1,2,4,8,16];
batchsizes = [1,2,4,8];
uniform_solvers = [1,2,3,4]; dynamic_solvers = 5 + uniform_solvers;
ns = length(uniform_solvers);  
test_function = 'rosenbrock';
solver = 'pounders';
titles = {'Balanced', 'Progressive', 'Imbalanced'};

width = 25;

legend_str = {};
for j = 1:ns
   legend_str{2*j-1} = strcat('Uniform, $r=$',num2str(batchsizes(j)));
   legend_str{2*j} = strcat('Dynamic, $r=$',num2str(batchsizes(j)));
end

figure; hlt = tiledlayout(1,3);
for j=1:3
    experiment_type = types{j};
    masterHf = zeros(num_probs,2*ns,num_epochs*m*(m+1));
    masterHg = zeros(num_probs,2*ns,num_epochs*m*(m+1));
    for macro_seed = 1:num_macro_seeds
        for micro_seed = 1:num_micro_seeds
            filename = strcat('results/',test_function,'_',experiment_type,'_',solver,'_',num2str(macro_seed),'_',num2str(micro_seed),'_',num2str(m),'.mat');
            load(filename);
            masterHf(num_micro_seeds*(macro_seed-1) + micro_seed,1:ns,:) = Hf(uniform_solvers,:);
            masterHf(num_micro_seeds*(macro_seed-1) + micro_seed,ns+1:2*ns,:) = Hf(dynamic_solvers,:);
            masterHg(num_micro_seeds*(macro_seed-1) + micro_seed,1:ns,:) = Hg(uniform_solvers,:);
            masterHg(num_micro_seeds*(macro_seed-1) + micro_seed,ns+1:2*ns,:) = Hg(dynamic_solvers,:);
        end
    end
    
    if strcmp(convergence,'fval')
        masterHf = permute(masterHf,[3 1 2]);
        prob_min = min(min(masterHf),[],3);
        H = masterHf;
    elseif strcmp(convergence,'gval')
        masterHg = permute(masterHg,[3 1 2]);
        prob_min = min(min(masterHg),[],3);
        H = masterHg; 
    end
    nexttile;
    [hlmu,hlmd] = plot_trajectory_percentiles2(H,prob_min,1:ns,ns+1:2*ns,width,mode,batchsizes,m);
    hold on;
    if j == 1
        hlm = [];
        for s = 1:ns
            hlm = cat(1,hlm,hlmu(s)); hlm = cat(1,hlm,hlmd(s));
        end
        legend(hlm,legend_str,'interpreter','latex','FontSize',FS);
    end
    if j == 1
        ylabel('$f(\mathbf{x}^k) - f(\mathbf{x}^*)$','interpreter','latex','FontSize',FS);
    end
    if j == 2
        if strcmp(mode,'eval')
            xlabel('Effective data passes','FontSize',FS);
        elseif strcmp(mode,'batch')
            xlabel('Batches of $r$ evaluations','FontSize',FS,'interpreter','latex');
        end
    end
    ylim([1e-16,1e5])
    xlim([0,1100])

    title(titles{j},'FontSize',FS);
end
hlt.Padding = 'compact'; hlt.TileSpacing = 'compact'; 
pause() % MAKE IT CAMERA READY DURING THE PAUSE
saveimgname = strcat('images/rosenbrock_uni_v_dyn_',test_function,'_','.eps');
saveas(gcf,saveimgname,'epsc');
end

%% EXPERIMENT 3 PLOTS
if ismember(3,which_plots)
% SPECIFIC TO EXPERIMENT 3: 
m = 8;
batchsizes = [1,2,4];
uniform_solvers = [1,2,3]; dynamic_solvers = 4 + uniform_solvers;
ns = length(uniform_solvers);  
test_function = 'cube';
solver = 'pounders';
titles = {'Balanced', 'Progressive', 'Imbalanced'};

width = 25;

legend_str = {};
for j = 1:ns
   legend_str{2*j-1} = strcat('Uniform, $r=$',num2str(batchsizes(j)));
   legend_str{2*j} = strcat('Dynamic, $r=$',num2str(batchsizes(j)));
end

figure; hlt = tiledlayout(1,3);
for j=1:3
    experiment_type = types{j};
    masterHf = NaN*ones(num_probs,2*ns,num_epochs*m*(m+1));
    masterHg = NaN*ones(num_probs,2*ns,num_epochs*m*(m+1));
    for macro_seed = 1:num_macro_seeds
        for micro_seed = 1:num_micro_seeds
            filename = strcat('results/',test_function,'_',experiment_type,'_',solver,'_',num2str(macro_seed),'_',num2str(micro_seed),'_',num2str(m),'.mat');
            try
                load(filename);
                masterHf(num_micro_seeds*(macro_seed-1) + micro_seed,1:ns,:) = Hf(uniform_solvers,:);
                masterHf(num_micro_seeds*(macro_seed-1) + micro_seed,ns+1:2*ns,:) = Hf(dynamic_solvers,:);
                masterHg(num_micro_seeds*(macro_seed-1) + micro_seed,1:ns,:) = Hg(uniform_solvers,:);
                masterHg(num_micro_seeds*(macro_seed-1) + micro_seed,ns+1:2*ns,:) = Hg(dynamic_solvers,:);
            catch
                fprintf(strcat(filename,' does not exist.\n'));
            end
        end
    end
    
    if strcmp(convergence,'fval')
        masterHf = permute(masterHf,[3 1 2]);
        prob_min = min(min(masterHf),[],3);
        H = masterHf;
    elseif strcmp(convergence,'gval')
        masterHg = permute(masterHg,[3 1 2]);
        prob_min = min(min(masterHg),[],3);
        H = masterHg; 
    end
    nexttile;
    [hlmu,hlmd] = plot_trajectory_percentiles2(H,prob_min,1:ns,ns+1:2*ns,width,mode,batchsizes,m);
    hold on;
    if j == 1
        hlm = [];
        for s = 1:ns
            hlm = cat(1,hlm,hlmu(s)); hlm = cat(1,hlm,hlmd(s));
        end
        legend(hlm,legend_str,'interpreter','latex','FontSize',FS);
    end
    if j == 1
        ylabel('$f(\mathbf{x}^k) - f(\mathbf{x}^*)$','interpreter','latex','FontSize',FS);
    end
    if j == 2
        if strcmp(mode,'eval')
            xlabel('Effective data passes','FontSize',FS);
        elseif strcmp(mode,'batch')
            xlabel('Batches of $r$ evaluations','FontSize',FS,'interpreter','latex');
        end
    end
    ylim([1e-16,1e5])
    xlim([0,900])
    title(titles{j},'FontSize',FS);
end
hlt.Padding = 'compact'; hlt.TileSpacing = 'compact'; 
pause() % MAKE IT CAMERA READY DURING THE PAUSE
saveimgname = strcat('images/cube_uni_v_dyn_',test_function,'_','.eps');
saveas(gcf,saveimgname,'epsc');
end

%% EXPERIMENT 4
if ismember(4,which_plots)
% SPECIFIC TO EXPERIMENT 1: 
m = 256;
batchsizes = [1,2,4,8,16,32,64,128,256];
num_batchsizes = length(batchsizes); 
test_function = 'logistic';
solver = 'fo';
titles = {'Balanced datasets', 'Progressive datasets', 'Imbalanced datasets'};
legendlabels = {};
for j = 1:num_batchsizes
    legendlabels{j} = strcat('SAM-FO, r=',num2str(batchsizes(j)));
end
num_gates = length(rup_gates); 
for k = 1:num_gates
    rup_gate = rup_gates(k);
    figure;
    hlt = tiledlayout(1,3);
    for j=1:3
        experiment_type = types{j};
       
        % RESOURCE UTILIZATION PLOT
        nexttile;
        [blh,medtable] = resource_utilization_plot3([num_batchsizes+1:2*num_batchsizes-1,num_batchsizes],num_macro_seeds,num_micro_seeds,num_epochs,m,test_function,experiment_type,solver,convergence,batchsizes,batchsizes,rup_gate);
        if j == 1
            ystr = '$\log_{2}\left(\displaystyle med_{\pi} R_{r,\pi,\tau,\mu}\right)$';
            ylabel(ystr,'interpreter','latex','FontSize',FS+2)
            legend(legendlabels,'Location','SouthWest');
        end
        if j == 2
           xlabel('$\mu$','interpreter','latex','FontSize',FS+2) 
        end
        title(titles{j});

        % CSV
        tablename = strcat('results/rup_',test_function,'_',experiment_type,'_',num2str(rup_gate),'.csv');
        writematrix(medtable,tablename);
    end
    hlt.Padding = 'compact'; hlt.TileSpacing = 'compact'; 
    title(hlt,strcat(rup_title,num2str(rup_gate)),'FontSize',FS,'interpreter','latex');
    pause() % MAKE IT CAMERA READY DURING THE PAUSE
    saveimgname = strcat('images/logistic_det_compare_',num2str(rup_gate),'.eps');
    saveas(gcf,saveimgname,'epsc');

end
end
%% EXPERIMENT 5
if ismember(5,which_plots)
% SPECIFIC TO EXPERIMENT 2: 
m = 16;
batchsizes = [1,2,4,8,16];
num_batchsizes = length(batchsizes); 
test_function = 'rosenbrock';
solver = 'pounders';
titles = {'Balanced', 'Progressive', 'Imbalanced'};
legendlabels = {};
for j = 1:num_batchsizes
    legendlabels{j} = strcat('SAM-POUNDERS, r=',num2str(batchsizes(j)));
end
num_gates = length(rup_gates); 
for k = 1:num_gates
    rup_gate = rup_gates(k);
    figure;
    hlt = tiledlayout(1,3);
    for j=1:3
        experiment_type = types{j};
        % RESOURCE UTILIZATION PLOT
        nexttile;
        [blh,medtable] = resource_utilization_plot3([num_batchsizes+1:2*num_batchsizes-1,num_batchsizes],num_macro_seeds,num_micro_seeds,num_epochs,m,test_function,experiment_type,solver,convergence,batchsizes,batchsizes,rup_gate);
        if j == 1
            ystr = '$\log_{2}\left(\displaystyle med_{\pi} R_{r,\pi,\tau,\mu}\right)$';
            ylabel(ystr,'interpreter','latex','FontSize',FS+2)
            legend(legendlabels,'Location','SouthWest');
        end
        if j == 2
           xlabel('$\mu$','interpreter','latex','FontSize',FS+2) 
        end
        title(titles{j});

        % CSV
        tablename = strcat('results/rup_',test_function,'_',experiment_type,'_',num2str(rup_gate),'.csv');
        writematrix(medtable,tablename);
    end
    hlt.Padding = 'compact'; hlt.TileSpacing = 'compact'; 
    title(hlt,strcat(rup_title,num2str(rup_gate)),'FontSize',FS,'interpreter','latex');
    pause() % MAKE IT CAMERA READY DURING THE PAUSE
    saveimgname = strcat('images/rosenbrock_det_compare_',num2str(rup_gate),'.eps');
    saveas(gcf,saveimgname,'epsc');
end
end

%% EXPERIMENT 6 PLOTS
if ismember(6,which_plots)
% SPECIFIC TO EXPERIMENT 3: 
m = 8;
batchsizes = [1,2,4,8];
num_batchsizes = length(batchsizes); 
test_function = 'cube';
solver = 'pounders';
titles = {'Balanced', 'Progressive', 'Imbalanced'};
legendlabels = {};
for j = 1:num_batchsizes
    legendlabels{j} = strcat('SAM-POUNDERS, r=',num2str(batchsizes(j)));
end
num_gates = length(rup_gates); 
for k = 1:num_gates
    rup_gate = rup_gates(k);
    figure;
    hlt = tiledlayout(1,3);
    for j=1:3
        experiment_type = types{j};
        % RESOURCE UTILIZATION PLOT
        nexttile;
        [blh,medtable] = resource_utilization_plot3([num_batchsizes+1:2*num_batchsizes-1,num_batchsizes],num_macro_seeds,num_micro_seeds,num_epochs,m,test_function,experiment_type,solver,convergence,batchsizes,batchsizes,rup_gate);
        if j == 1
            ystr = '$\log_{2}\left(\displaystyle med_{\pi} R_{r,\pi,\tau,\mu}\right)$';
            ylabel(ystr,'interpreter','latex','FontSize',FS+2)
            legend(legendlabels,'Location','SouthWest');
        end
        if j == 2
           xlabel('$\mu$','interpreter','latex','FontSize',FS+2) 
        end
        title(titles{j});

        % CSV
        tablename = strcat('results/rup_',test_function,'_',experiment_type,'_',num2str(rup_gate),'.csv');
        writematrix(medtable,tablename);
    end
    hlt.Padding = 'compact'; hlt.TileSpacing = 'compact'; 
    title(hlt,strcat(rup_title,num2str(rup_gate)),'FontSize',FS,'interpreter','latex');
    pause() % MAKE IT CAMERA READY DURING THE PAUSE
    saveimgname = strcat('images/cube_det_compare_',num2str(rup_gate),'.eps');
    saveas(gcf,saveimgname,'epsc');
end
end

%% FIGURE 7 PLOTS
if ismember(7,which_plots)
% SPECIFIC TO EXPERIMENT 1: 
m = 256;
batchsizes = [1,2,4,8,16,32,64,128,256];
num_batchsizes = length(batchsizes); 
test_function = 'logistic';
solver = 'fo';
titles = {'Balanced datasets', 'Progressive datasets', 'Imbalanced datasets'};

for j=1:3
    experiment_type = types{j};
    masterHf = zeros(num_probs,2*(num_batchsizes-1),num_epochs*m);
    masterHg = zeros(num_probs,2*(num_batchsizes-1),num_epochs*m);
    for macro_seed = 1:num_macro_seeds
        for micro_seed = 1:num_micro_seeds
            filename = strcat('results/',test_function,'_',experiment_type,'_',solver,'_',num2str(macro_seed),'_',num2str(micro_seed),'_',num2str(m),'.mat');
            load(filename);
            masterHf(num_micro_seeds*(macro_seed-1) + micro_seed,1:num_batchsizes-1,:) = Hf(num_batchsizes+1:2*num_batchsizes-1,:);
            masterHg(num_micro_seeds*(macro_seed-1) + micro_seed,1:num_batchsizes-1,:) = Hg(num_batchsizes+1:2*num_batchsizes-1,:);
            filename = strcat('results/no_lip_',test_function,'_',experiment_type,'_',solver,'_',num2str(macro_seed),'_',num2str(micro_seed),'_',num2str(m),'.mat');
            load(filename);
            masterHf(num_micro_seeds*(macro_seed-1) + micro_seed,num_batchsizes:2*(num_batchsizes-1),:) = Hf(num_batchsizes+1:2*num_batchsizes-1,:);
            masterHg(num_micro_seeds*(macro_seed-1) + micro_seed,num_batchsizes:2*(num_batchsizes-1),:) = Hg(num_batchsizes+1:2*num_batchsizes-1,:);
        end
    end
    
    if strcmp(convergence,'fval')
        masterHf = permute(masterHf,[3 1 2]);
        prob_min = min(min(masterHf),[],3);
        H = masterHf;
    elseif strcmp(convergence,'gval')
        masterHg = permute(masterHg,[3 1 2]);
        prob_min = min(min(masterHg),[],3);
        H = masterHg; 
    end
    figure;
    hlt = tiledlayout(1,4);
    for b = 1:2:(num_batchsizes-1)
        nexttile
        [hl,~] = plot_trajectory_percentiles(H,prob_min,[b,num_batchsizes-1+b],m);
        hold on
        %title(strcat(test_function,', ',experiment_type,' ,batchsize=',num2str(batchsizes(b))));
        if b == 1
            legend(hl, {'SAM-FO with Lipschitz constants','SAM-FO without Lipschitz constants'},'FontSize',FS)
        end
        if b == 1
            ylabel('$f(\mathbf{x}^k) - f(\mathbf{x}^*)$','interpreter','latex','FontSize',FS);
        end
        if b == 1
            xlabel('Effective data passes','FontSize',FS);
        end
        if j == 1
            ylim([1e-17,1e5])
        elseif j == 2
            ylim([1e-10,1e2])
        elseif j == 3
            ylim([1e-16,1e2])
        end
        title(strcat('r =  ',num2str(batchsizes(b))),'FontSize',FS);
    end
    hlt.Padding = 'compact'; hlt.TileSpacing = 'compact'; 
    title(hlt,titles{j},'FontSize',FS);
    pause() % MAKE IT CAMERA READY DURING THE PAUSE
    saveimgname = strcat('images/lip_v_nolip_',test_function,'_',experiment_type,'.eps');
    saveas(gcf,saveimgname,'epsc');
end

    
end

%% FIGURE 8 PLOTS
if ismember(8,which_plots)
% SPECIFIC TO EXPERIMENT 2: 
m = 16;
batchsizes = [1,2,4,8,16];
num_batchsizes = length(batchsizes); 
test_function = 'rosenbrock';
solver = 'pounders';
titles = {'Balanced', 'Progressive', 'Imbalanced'};

for j=1:3
    experiment_type = types{j};
    masterHf = zeros(num_probs,2*(num_batchsizes-1),num_epochs*m*(m+1));
    masterHg = zeros(num_probs,2*(num_batchsizes-1),num_epochs*m*(m+1));
    for macro_seed = 1:num_macro_seeds
        for micro_seed = 1:num_micro_seeds
            filename = strcat('results/',test_function,'_',experiment_type,'_',solver,'_',num2str(macro_seed),'_',num2str(micro_seed),'_',num2str(m),'.mat');
            load(filename);
            masterHf(num_micro_seeds*(macro_seed-1) + micro_seed,1:num_batchsizes-1,:) = Hf(num_batchsizes+1:2*num_batchsizes-1,:);
            masterHg(num_micro_seeds*(macro_seed-1) + micro_seed,1:num_batchsizes-1,:) = Hg(num_batchsizes+1:2*num_batchsizes-1,:);
            filename = strcat('results/no_lip_',test_function,'_',experiment_type,'_',solver,'_',num2str(macro_seed),'_',num2str(micro_seed),'_',num2str(m),'.mat');
            load(filename);
            masterHf(num_micro_seeds*(macro_seed-1) + micro_seed,num_batchsizes:2*(num_batchsizes-1),:) = Hf(num_batchsizes+1:2*num_batchsizes-1,:);
            masterHg(num_micro_seeds*(macro_seed-1) + micro_seed,num_batchsizes:2*(num_batchsizes-1),:) = Hg(num_batchsizes+1:2*num_batchsizes-1,:);
        end
    end
    
    if strcmp(convergence,'fval')
        masterHf = permute(masterHf,[3 1 2]);
        prob_min = min(min(masterHf),[],3);
        H = masterHf;
    elseif strcmp(convergence,'gval')
        masterHg = permute(masterHg,[3 1 2]);
        prob_min = min(min(masterHg),[],3);
        H = masterHg; 
    end
    figure;
    hlt = tiledlayout(1,4);
    for b = 1:(num_batchsizes-1)
        nexttile
        [hl,~] = plot_trajectory_percentiles(H,prob_min,[b,num_batchsizes-1+b],m);
        hold on
        %title(strcat(test_function,', ',experiment_type,' ,batchsize=',num2str(batchsizes(b))));
        if b == 1
            legend(hl, {'SAM-POUNDERS with Lipschitz constants','SAM-POUNDERS without Lipschitz constants'},'FontSize',FS)
        end
        if b == 1
            ylabel('$f(\mathbf{x}^k) - f(\mathbf{x}^*)$','interpreter','latex','FontSize',FS);
        end
        if b == 1
            xlabel('Effective data passes','FontSize',FS);
        end
        ylim([1e-16,1e5])
        xlim([0,1100])
        title(strcat('r=',num2str(batchsizes(b))),'FontSize',FS);
    end
    hlt.Padding = 'compact'; hlt.TileSpacing = 'compact'; 
    title(hlt,titles{j},'FontSize',FS);
    pause() % MAKE IT CAMERA READY DURING THE PAUSE
    saveimgname = strcat('images/lip_v_nolip_',test_function,'_',experiment_type,'.eps');
    saveas(gcf,saveimgname,'epsc');
end

end

%% FIGURE 9 PLOTS
if ismember(9,which_plots)
% SPECIFIC TO EXPERIMENT 3: 
m = 8;
batchsizes = [1,2,4,8];
num_batchsizes = length(batchsizes); 
test_function = 'cube';
solver = 'pounders';
titles = {'Balanced', 'Progressive', 'Imbalanced'};

for j=1:3
    experiment_type = types{j};
    masterHf = NaN*ones(num_probs,2*(num_batchsizes-1),num_epochs*m*(m+1));
    masterHg = NaN*ones(num_probs,2*(num_batchsizes-1),num_epochs*m*(m+1));
    for macro_seed = 1:num_macro_seeds
        for micro_seed = 1:num_micro_seeds
            filename = strcat('results/',test_function,'_',experiment_type,'_',solver,'_',num2str(macro_seed),'_',num2str(micro_seed),'_',num2str(m),'.mat');
            load(filename);
            masterHf(num_micro_seeds*(macro_seed-1) + micro_seed,1:num_batchsizes-1,:) = Hf(num_batchsizes+1:2*num_batchsizes-1,:);
            masterHg(num_micro_seeds*(macro_seed-1) + micro_seed,1:num_batchsizes-1,:) = Hg(num_batchsizes+1:2*num_batchsizes-1,:);
            filename = strcat('results/no_lip_',test_function,'_',experiment_type,'_',solver,'_',num2str(macro_seed),'_',num2str(micro_seed),'_',num2str(m),'.mat');
            load(filename);
            masterHf(num_micro_seeds*(macro_seed-1) + micro_seed,num_batchsizes:2*(num_batchsizes-1),:) = Hf(num_batchsizes+1:2*num_batchsizes-1,:);
            masterHg(num_micro_seeds*(macro_seed-1) + micro_seed,num_batchsizes:2*(num_batchsizes-1),:) = Hg(num_batchsizes+1:2*num_batchsizes-1,:);
        end
    end
    
    if strcmp(convergence,'fval')
        masterHf = permute(masterHf,[3 1 2]);
        prob_min = min(min(masterHf),[],3);
        H = masterHf;
    elseif strcmp(convergence,'gval')
        masterHg = permute(masterHg,[3 1 2]);
        prob_min = min(min(masterHg),[],3);
        H = masterHg; 
    end
    figure;
    hlt = tiledlayout(1,3);
    for b = 1:(num_batchsizes-1)
        nexttile
        [hl,~] = plot_trajectory_percentiles(H,prob_min,[b,num_batchsizes-1+b],m);
        hold on
        %title(strcat(test_function,', ',experiment_type,' ,batchsize=',num2str(batchsizes(b))));
        if b == 1
            legend(hl, {'SAM-POUNDERS with Lipschitz constants','SAM-POUNDERS without Lipschitz constants'},'FontSize',FS)
        end
        if b == 1
            ylabel('$f(\mathbf{x}^k) - f(\mathbf{x}^*)$','interpreter','latex','FontSize',FS);
        end
        if b == 2
            xlabel('Effective data passes','FontSize',FS);
        end
        ylim([1e-16,1e5])
        xlim([0,800])
        title(strcat('r=',num2str(batchsizes(b))),'FontSize',FS);
    end
    hlt.Padding = 'compact'; hlt.TileSpacing = 'compact'; 
    title(hlt,titles{j},'FontSize',FS);
    pause() % MAKE IT CAMERA READY DURING THE PAUSE
    saveimgname = strcat('images/lip_v_nolip_',test_function,'_',experiment_type,'.eps');
    saveas(gcf,saveimgname,'epsc');
end

end

%% FIGURE 10 PLOTS
if ismember(10,which_plots)
% SPECIFIC TO EXPERIMENT 10: 
m = 256;
num_epochs = 100;

titles = {'Balanced datasets', 'Progressive datasets', 'Imbalanced datasets'};

figure; hlt = tiledlayout(1,3);
for j=1:3
    experiment_type = types{j};
    masterHf = zeros(num_probs,2,num_epochs*m);
    masterHg = zeros(num_probs,2,num_epochs*m);
    for macro_seed = 1:num_macro_seeds
        for micro_seed = 1:num_micro_seeds
            filename = strcat('results/sag_',experiment_type,'_',num2str(macro_seed),'_',num2str(micro_seed),'_',num2str(m),'.mat');
            load(filename);
            masterHf(num_micro_seeds*(macro_seed-1) + micro_seed,:,:) = Hf;
            masterHg(num_micro_seeds*(macro_seed-1) + micro_seed,:,:) = Hg;
        end
    end
    
    if strcmp(convergence,'fval')
        masterHf = permute(masterHf,[3 1 2]);
        prob_min = min(min(masterHf),[],3);
        H = masterHf;
    elseif strcmp(convergence,'gval')
        masterHg = permute(masterHg,[3 1 2]);
        prob_min = min(min(masterHg),[],3);
        H = masterHg; 
    end
    
    nexttile;
    [hl,~] = plot_trajectory_percentiles(H,prob_min,[2 1],m);
    hold on
    title(titles{j},'FontSize',FS);
    if j == 1
        legend(hl,{'SAG-LS (Lipschitz)','Dynamic SAM-FO'},'FontSize',FS);
        ylabel('$f(\mathbf{x}^k) - f(\mathbf{x}^*)$','interpreter','latex','FontSize',FS);
    end
    if j == 2
        xlabel('Effective data passes','FontSize',FS);
    end
end
    pause(); % MAKE CAMERA READY
    saveimgname = strcat('images/sag.eps');
    saveas(gcf,saveimgname,'epsc');
end

if ismember(11,which_plots)   
    b = 1; micro_seed = 1; macro_seed = 1; num_epochs = 100; m = 256;

    titles = {'Balanced generation', 'Progressive generation', 'Imbalanced generation'};
    types = {'balanced','progressive','imbalanced'};
    test_function = 'logistic';
    for j = 1:3
        hlt = just_one_run(test_function,types{j},b,micro_seed,macro_seed,num_epochs,m);
        title(hlt,titles{j},'FontSize',FS);   
        hlt.Padding = 'compact'; hlt.TileSpacing = 'compact'; 
        ylim([0.5 m+0.5]);
        pause() %MAKE CAMERA READY
        saveimgname = strcat('images/',types{j},'_',test_function,'_visualize.eps');
        saveas(gcf,saveimgname,'epsc');
    end

    
end

if ismember(12,which_plots)   
    b = 1; micro_seed = 1; macro_seed = 1; num_epochs = 100; m = 16;

    titles = {'Balanced generation', 'Progressive generation', 'Imbalanced generation'};
    types = {'balanced','progressive','imbalanced'};
    test_function = 'rosenbrock';
    for j = 1:3
        hlt = just_one_run(test_function,types{j},b,micro_seed,macro_seed,num_epochs,m);
        title(hlt,titles{j},'FontSize',FS);   
        hlt.Padding = 'compact'; hlt.TileSpacing = 'compact'; 
        ylim([0.5 m+0.5]);
        pause() %MAKE CAMERA READY
        saveimgname = strcat('images/',types{j},'_',test_function,'_visualize.eps');
        saveas(gcf,saveimgname,'epsc');
    end

    
end

if ismember(13,which_plots)   
    b = 1; micro_seed = 1; macro_seed = 1; num_epochs = 100; m = 8;

    titles = {'Balanced generation', 'Progressive generation', 'Imbalanced generation'};
    types = {'balanced','progressive','imbalanced'};
    test_function = 'cube';
    for j = 1:3
        hlt = just_one_run(test_function,types{j},b,micro_seed,macro_seed,num_epochs,m);
        title(hlt,titles{j},'FontSize',FS);   
        hlt.Padding = 'compact'; hlt.TileSpacing = 'compact'; 
        ylim([0.5 m+0.5]);
        pause() %MAKE CAMERA READY
        saveimgname = strcat('images/',types{j},'_',test_function,'_visualize.eps');
        saveas(gcf,saveimgname,'epsc');
    end

    
end
end