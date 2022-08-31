function [Lip,hlt] = just_one_run(test_function,data_type,macro_seed,micro_seed,b,num_epochs,m)

addpath('~/sampounders/');
global spsolver

nonzerocolor = [0 0.447 0.741];
logcolor = [0.85 0.325 0.098];

FS = 16; LW = 3;

n = m;
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
    alpha((m/2)+1:m) = ones(m/2,1); %alpha_saved(2:2:m);
    lipY = zeros(1,m); 
    lipY(1:m/2) = (alpha(1:m/2))*20;
    func = @(x,Set)generalized_rosenbrock(x,Set,alpha);
    solver = 'pounders';
elseif strcmp(test_function,'cube')
    C = 1;
    X0 = 2*(-0.5 + rand(1,n));
    lipY = 30*C*alpha; lipY(1) = 0;
    func = @(x,Set)generalized_cube(x,Set,alpha);
    solver = 'pounders';
elseif strcmp(test_function,'logistic')
    addpath(genpath('SGDLibrary-master/')); 
    alpha = 4*alpha/(m);
    data = logistic_regression_data_generator(m, n, alpha);
    X0 = data.w_init';
    lambda = 1e-1; % just to make strongly convex
    lipY = (sum(0.25*(data.x_train).^2) + lambda)/m;
    %lipY = lambda*ones(m,1)/m; 
    func = @(x,Set)log_reg(x,Set,data.x_train,data.y_train,lambda);
    solver = 'fo';
else
    error('test_function must be set to rosenbrock, cube, or logistic.');
end

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
spsolver=2; addpath('../minq5/');

figure;
hlt = tiledlayout(1,2); nexttile;
yyaxis left
%ax1 = axes(hlt, 'NextPlot', 'add', 'YAxisLocation', 'right', 'Box', 'on');
rng(micro_seed); % Controls randomness of run (as opposed to generation of problem data) 
if strcmp(solver,'fo')
    [X,~,~,~,Eval,Lip] = sam_fo(func,X0,n,npmax,nfmax,gtol,delta,nfs,m,F0,xind,Low,Upp,printf,b,'adaptive',lipY);
    a = sum(Eval,2);
    iters = nnz(a);
    Lip = Lip(1:iters,:);
    hlspy = imshow(~Eval(1:iters,:)',[nonzerocolor;1.0 1.0 1.0])%,'Parent',ax1);
    axis on;
    totals = sum(Eval,1);
    
elseif strcmp(solver,'pounders')
    [X,~,~,~,Eval,Lip] = sam_pounders(func,X0,n,npmax,nfmax,gtol,delta,nfs,m,F0,xind,Low,Upp,printf,b,'adaptive');%,lipY);
    if strcmp(test_function,'rosenbrock')
        % shuffle to match the paper
        savedEval = Eval;
        Eval(:,2:2:m) = savedEval(:,((m/2)+1):m);
        Eval(:,1:2:(m-1)) = savedEval(:,1:(m/2));
    end
    a = sum(Eval,2);
    iters = nnz(a);
    %spacing = 1:m:iters;
    %Eval2 = zeros(m,length(spacing));
    %for j = 1:length(spacing)
    %    Eval2(:,j) = sum(Eval(1+(j-1)*m:j*m,:),1)';
    %end
    Lip = Lip(1:iters,:);
    hlspy = imshow(~Eval(1:iters,:)',[nonzerocolor;1.0 1.0 1.0])%,'Parent',ax1);
    %hlspy = imshow(~Eval2,[nonzerocolor;1.0 1.0 1.0]);
    axis on;
    totals = sum(Eval,1); 
end

ax1 = gca;
ax1.YColor = 'black';
ylabel('Component function number','FontSize',FS);

fvals = zeros(1,iters);
for j = 1:iters
    try
        F = func(X(j,:),1:m);
        if strcmp(solver,'fo')
            fvals(j) = sum(F);
        elseif strcmp(solver,'pounders')
            fvals(j) = min(sum(F.^2));
        end
    catch
        % there was an indexing issue
        fvals(end) = [];
        iters = iters - 1;
    end
end

if strcmp(test_function,'logistic')
    [fmin,imin] = min(fvals);
    replace_fmin = fvals(imin-1);
else
    fmin = 0;
    replace_fmin = min(fvals);
end

% smooth
for j = 2:iters
    if fvals(j) == fmin
        fvals(j) = replace_fmin;
    end
    fvals(j) = min(fvals(j),fvals(j-1));
end
    
yyaxis right
%ax2 = axes(hlt, 'NextPlot', 'add');
hl = semilogy(fvals-fmin,'LineWidth',LW,'LineStyle','--','Color',logcolor);
ax2 = gca;
xlim(ax2,[1 iters]);
ax2.YColor = 'black'; 
%ylabel('$f(\mathbf{x}^k) - f(\mathbf{x}^*)$','interpreter','latex','FontSize',FS);
ax2.XColor = 'black';
xlabel('k','Interpreter','latex','FontSize',FS); 

nexttile;
barh(1:m,totals,'FaceColor',nonzerocolor);
set(gca,'YDir','reverse')
if strcmp(solver,'fo')
    xlabel('Gradient evaluations of component','FontSize',FS);
else
    xlabel('Function evaluations of component','FontSize',FS);
end

end