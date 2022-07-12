function experiment3(data_type, macro_seed, micro_seed, num_epochs, m)

% RUNS SAG COMPARISONS. 

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

try
    addpath(genpath('SAG/')); 
catch
    error('Please install SAG.');
end
alpha = 4*alpha/(m);
data = logistic_regression_data_generator(m, n, alpha);
X0 = data.w_init';
lambda = 1e-1; % just to make strongly convex
lipY = (sum(0.25*(data.x_train).^2) + lambda)/m;
func = @(x,Set)log_reg(x,Set,data.x_train,data.y_train,lambda);

%% FIXED THINGS FOR SOLVERS
%npmax = 2*n+1;
npmax = (n+1)*(n+2)/2;
nfmax = num_epochs*m;
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
Hf = NaN*ones(2,nfmax);
Hg = NaN*ones(2,nfmax);

%% First run the solver using adaptive batchsize generation of I_k, J_k

% only testing batchsize 1: 
b = 1; 

rng(micro_seed); % Controls randomness of run (as opposed to generation of problem data) 
[X,~,~,xkin,Eval] = sam_fo(func,X0,n,npmax,nfmax,gtol,delta,nfs,m,F0,xind,Low,Upp,printf,b,'adaptive',lipY);

% Post processing
gradcount = cumsum(sum(Eval(1:xkin,:),2))';
fvals = zeros(1,xkin); gvals = zeros(1,xkin);
for j = 1:xkin
    [F,J] = func(X(j,:),1:m);
    fvals(j) = sum(F);
    gvals(j) = norm(sum(J));
end

% Save to output data
prev = 0;
for kk = 1:length(gradcount)
    if gradcount(kk) > nfmax
        break
    end
    Hf(1,prev + 1:gradcount(kk)) = fvals(kk)*ones(1,gradcount(kk)-prev);
    Hg(1,prev + 1:gradcount(kk)) = gvals(kk)*ones(1,gradcount(kk)-prev);
    prev = gradcount(kk); 
end


%% Now run SAG

rng(micro_seed);
% specific things
X = data.x_train; y = data.y_train;

maxIter = m*num_epochs;

Xt = X'; % Function works with transpose of X

xtx = sum(Xt.^2,2);

Lmax = max(lipY);
randVals = rand(maxIter,2); % Random values to be used by the algorithm

d = zeros(n,1);
g = zeros(m,1);
covered = int32(zeros(m,1));

history = zeros(1,maxIter*m);
SAG_LipschitzLS_logistic(X0',X,y',lambda,Lmax,m*lipY',randVals,d,g,covered,int32(1),xtx,history);

% Post processing
fvals = zeros(1,maxIter); gvals = zeros(1,maxIter);
for j = 1:maxIter
    [F,J] = func(history(n*(j-1)+1:n*j),1:m);
    fvals(j) = sum(F);
    gvals(j) = norm(sum(J));
end

% Save to output data
for kk = 1:num_epochs
    Hf(2,((kk-1)*m + 1):kk*m) = fvals(kk)*ones(1,m);
    Hg(2,((kk-1)*m + 1):kk*m) = gvals(kk)*ones(1,m);
end


% Save output data to file in results directory.
filename = strcat('results/sag_',data_type,'_',num2str(macro_seed),'_',num2str(micro_seed),'_',num2str(m),'.mat');
save(filename,'Hf','Hg'); 

end
