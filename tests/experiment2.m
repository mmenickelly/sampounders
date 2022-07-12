function experiment2(test_function, data_type, solver, macro_seed, micro_seed, batchsizes, num_epochs, m)


% RUNS LIPSCHITZ VS NO LIPSCHITZ TESTS. 

global spsolver

n = m;

%% Generate one random instance of the problem with fixed seed.
rng(macro_seed); 
if strcmp(data_type,'imbalanced') 
   alpha =  ones(1,m); alpha(m) = m; alpha(m-1) = m;
elseif strcmp(data_type,'progressive')
    alpha = 1:m;
elseif strcmp(data_type,'balanced')
    alpha = ones(1,m);
else
    error('data_type must be set to imbalanced, progressive, or balanced');
end

if strcmp(test_function,'rosenbrock')
    C = 1;
    X0 = 2*C*(-0.5 + rand(1,n));
    alpha_saved = alpha;
    alpha(1:m/2) = alpha_saved(1:2:m-1);
    alpha((m/2)+1:m) = alpha_saved(2:2:m);
    lipY = zeros(1,m); 
    lipY(1:m/2) = (alpha(1:m/2))*20;
    func = @(x,Set)generalized_rosenbrock(x,Set,alpha);
    if strcmp(solver,'fo')
        func = @(x,Set)square_each_component(func,x,Set);
    end
elseif strcmp(test_function,'cube')
    C = 1;
    X0 = 2*C*(-0.5 + rand(1,n));
    lipY = 3*C*alpha; lipY(1) = 0;
    func = @(x,Set)generalized_cube(x,Set,alpha);
    if strcmp(solver,'fo')
        func = @(x,Set)square_each_component(func,x,Set);
    end
elseif strcmp(test_function,'logistic')
    alpha = 4*alpha/(m);
    data = logistic_regression_data_generator(m, n, alpha);
    X0 = data.w_init';
    lambda = 1e-1; % just to make strongly convex
    lipY = (sum(0.25*(data.x_train).^2) + lambda)/m;
    func = @(x,Set)log_reg(x,Set,data.x_train,data.y_train,lambda);
else
    error('test_function must be set to rosenbrock, cube, or logistic.');
end

%% FIXED THINGS FOR SOLVERS
%npmax = 2*n+1;
npmax = (n+1)*(n+2)/2;
if strcmp(solver,'fo')
	nfmax = num_epochs*m;
else
	nfmax = num_epochs*m*(n+1);
end
gtol = 1e-13;
delta = 1;
nfs = 0;
F0 = [];
xind = 1;
Low = -Inf(1,n);
Upp = Inf(1,n);
printf = true;
%printf = false;

% Choose your solver: - should eventually be phased out in POptUS? 
spsolver=2; addpath('../minq5/');

%% Instantiate output data
num_batchsizes = length(batchsizes);
Hf = NaN*ones(2*num_batchsizes-1,nfmax);
Hg = NaN*ones(2*num_batchsizes-1,nfmax);
Evals = cell(1,2*num_batchsizes-1);

%% Run the solver using dynamic generation of I_k, J_k
for k = 1:(num_batchsizes - 1)
    b = batchsizes(k);
    rng(micro_seed);
    if strcmp(solver,'fo')
        [X,~,~,xkin,Eval] = sam_fo(func,X0,n,npmax,nfmax,gtol,delta,nfs,m,F0,xind,Low,Upp,printf,b,'adaptive');
    else
        [X,~,~,xkin,Eval] = sam_pounders(func,X0,n,npmax,nfmax,gtol,delta,nfs,m,F0,xind,Low,Upp,printf,b,'adaptive');
    end

    % Post processing
    gradcount = cumsum(sum(Eval(1:xkin,:),2))';
    fvals = zeros(1,xkin); gvals = zeros(1,xkin);
    for j = 1:xkin
        [F,J] = func(X(j,:),1:m);
        if strcmp(test_function,'logistic')
            fvals(j) = sum(F);
            gvals(j) = norm(sum(J));
        else
            fvals(j) = sum(F.^2);
            gvals(j) = norm(F*J);
        end
    end
    
    % Save to output data
    prev = 0;
    for kk = 1:length(gradcount)
        if gradcount(kk) > nfmax
            break
        end
        Hf(num_batchsizes + k,prev + 1:gradcount(kk)) = fvals(kk)*ones(1,gradcount(kk)-prev);
        Hg(num_batchsizes + k,prev + 1:gradcount(kk)) = gvals(kk)*ones(1,gradcount(kk)-prev);
        prev = gradcount(kk); 
    end
    Evals{num_batchsizes + k} = Eval;
end

% Save output data to file in results directory.
filename = strcat('results/no_lip_',test_function,'_',data_type,'_',solver,'_',num2str(macro_seed),'_',num2str(micro_seed),'_',num2str(m),'.mat');
save(filename,'Hf','Hg','Evals'); 

end
