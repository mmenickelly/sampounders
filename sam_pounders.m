function [X,F,flag,xkin,Eval,Lip] = ...
    sam_pounders(fun,X0,n,npmax,nfmax,gtol,delta,nfs,m,F0,xkin,L,U,printf,batchsize,sampling,lipY)

% Choose your solver:
global spsolver % Value is fixed outside
spsolver=2; addpath('minq5/'); % Arnold Neumaier's minq

% 0. Check inputs
[flag,X0,npmax,F0,L,U] = ...
    checkinputss(fun,X0,n,npmax,nfmax,gtol,delta,nfs,m,F0,xkin,L,U);
if (flag == -1), X=[]; F=[]; return; end % Problem with the input

% --INTERNAL PARAMETERS [won't be changed elsewhere, defaults in ( ) ]-----
maxdelta = min(.5*min(U-L),1e3*delta); % [dbl] Maximum tr radius
mindelta = min(delta*1e-13,gtol/10); % [dbl] Min tr radius (technically 0)
gam0 = .5;      % [dbl] Parameter in (0,1) for shrinking delta  (.5)
gam1 = 2;       % [dbl] Parameter >1 for enlarging delta   (2)
eta1 = .1;     % [dbl] Parameter 2 for accepting point, 0<eta1<1 (.2)
Par(1) = sqrt(n); % [dbl] delta multiplier for checking validity
Par(2) = max(10,sqrt(n)); % [dbl] delta multiplier for all interp. points
Par(3) = 1e-3;  % [dbl] Pivot threshold for validity (1e-5)
Par(4) = .001;  % [dbl] Pivot threshold for additional points (.001)
eta2 = 0; % [dbl] Parameter > 0 that determines acceptability
if printf
    if strcmp(sampling,'arbitrary') || strcmp(sampling,'adaptive')
        disp('  nf   delta           f0           g0           var_m           var_f');
        progstr = '%4i %9.2e  %11.5e %12.4e %11.5e %11.5e \n'; % Line-by-line
    else
        disp('  nf   delta           f0           g0           ');
        progstr = '%4i %9.2e  %11.5e %12.4e \n'; % Line-by-line
    end
end
% -------------------------------------------------------------------------

% --INTERMEDIATE VARIABLES-------------------------------------------------
% D       [dbl] [1-by-n] Generic displacement vector
% G       [dbl] [n-by-1] Model gradient at X(xkin,:)
% H       [dbl] [n-by-n] Model Hessian at X(xkin,:)
% Hdel    [dbl] [n-by-n] Change to model Hessian at X(xkin,:)
% Lows    [dbl] [1-by-n] Vector of subproblem lower bounds
% Upps    [dbl] [1-by-n] Vector of subproblem upper bounds
% Mdir    [dbl] [n-by-n] Unit row directions to improve model/geometry
% Mind    [int] [npmax-by-1] Integer vector of model interpolation indices
% Xsp     [dbl] [1-by-n] Subproblem solution
% c       [dbl] Model value at X(xkin,:)
% mdec    [dbl] Change predicted by the model, m(nf)-m(xkin)
% nf      [int] Counter for the number of function evaluations
% ng      [dbl] Norm of (projection of) G
% np      [int] Number of model interpolation points
% rho     [dbl] Ratio of actual decrease to model decrease
% valid   [log] Flag saying if model is fully linear within Par(1)*delta
% -------------------------------------------------------------------------

if nfs==0 % Need to do the first evaluation
    X = [X0; zeros(nfmax-1,n)]; % Stores the point locations
    F = nan(nfmax,m); % Stores the function values
    Eval = logical(zeros(nfmax,m)); % Stores whether eval i includes component j
    nf = 1;
    F(nf,:) = fun(X(nf,:),1:m);
    Eval(nf,1:m) = 1; 
    if printf
        fprintf('%4i    Initial point  %11.5e\n',nf,sum(F(nf,:).^2));
    end
else % Have other function values around
    X = [X0(1:max(1,nfs),:); zeros(nfmax,n)]; % Stores the point locations
    F = [F0(1:nfs,:); nan(nfmax,m)]; % Stores the function values
    %!! For now assume that all m components evaluated with incoming points
    Eval = logical([ones(nfs,m); zeros(nfmax,m)]);
    nf = nfs;
    nfmax = nfmax+nfs;
end

% For studying behavior of Lipschitz constant:
Lip = F;

Res = zeros(size(F)); % Stores the residuals for model updates
Fy = F(xkin,:); Gres=zeros(n,m); Hres = zeros(n,n,m);
center_ind = ones(1,m); % Stores model centers
old_deltas = delta*ones(1,m);

if nargin == 16
    lipY = ones(1,m); % default lipschitz constant estimates
    nolip = true; 
    first_success_now = false; 
    first_success_already = false;
else
    nolip = false;
    first_success_now = true;
    first_success_already = true;
end

Lip(1,:) = lipY; 

% evaluate all components immediately:
crit_check = true;

while sum(sum(Eval))<nfmax && delta > mindelta
    %% Step 1: Choose a subset
    if nargin == 16 && first_success_now && ~first_success_already
        crit_check = true;
        first_success_already = true; 
    end
    if crit_check
        subset_m = 1:m; probs_m = ones(1,m); var_m = 0;
        subset_f = 1:m; probs_f = ones(1,m); var_f = 0;
    else
        if strcmp(sampling,'arbitrary')
            [subset_m,probs_m,var_m] = arbitrary_selection_tr(X,Fy,lipY,batchsize,delta,old_deltas,xkin,center_ind);
        elseif strcmp(sampling,'adaptive')
            [subset_m,probs_m,var_m] = adaptive_selection_tr(X,Fy,lipY,batchsize,delta,old_deltas,xkin,center_ind);
        elseif strcmp(sampling,'uniform')
            [subset_m,probs_m] = uniform_selection(m,batchsize,xkin,center_ind,'center');
        else
            error('sampling must be arbitrary or uniform');
        end
    end
    
    %% Step 2: Update model centers and models
    % update center_inds at subset_model and evaluate those centers
    old_center_ind = center_ind;
    center_ind(subset_m) = xkin;
    F(xkin,subset_m) = fun(X(xkin,:),subset_m);
    
    Fy(subset_m) = F(xkin,subset_m);
    Eval(xkin,subset_m) = true;
    old_deltas(subset_m) = delta;
    
    [nf,Eval,Res,X,F,Hres,lipY,G,H,Gres,Lip] = ... 
        update_model_centers_and_models(fun,xkin,n,m,subset_m,delta,nf,Eval,center_ind,old_center_ind,Res,X,F,Fy,Gres,Hres,npmax,Par,L,U,nfmax,lipY,probs_m,printf,nolip,Lip);
        
    %% Step 2.5: Check for small gradient norm and acceptability
    
    ng = norm(G.*( and(X(xkin,:)>L,G'>0) + and(X(xkin,:)<U,G'<0) )');
    
    if ng < gtol && ~crit_check
        crit_check = true;
    elseif ng < gtol && crit_check
        disp('g is sufficiently small');
        flag = 0; return;
    else
        crit_check = false;
    end
    
    % Check for acceptability
    if ng <= eta2*sum(lipY)*delta 
        acceptable = false;
    else
        acceptable = true;
    end
    
    if acceptable
       %% Step 3: Solve TRSP to compute step
        Lows = max((L-X(xkin,:)),-delta);
        Upps = min((U-X(xkin,:)),delta);
        if spsolver==2 % Arnold Neumaier's minq
            [Xsp,mdec] = minqsw(0,G,H,Lows',Upps',0,zeros(n,1));
        end
        Xsp = Xsp'; % Solvers currently work with column vectors
        
       %% Step 4: Choose a second subset
        if strcmp(sampling,'arbitrary')
            [subset_f,probs_f,var_f] = arbitrary_selection_twopt(X,Xsp,Fy,lipY,batchsize,delta,old_deltas,xkin,center_ind,'trial');
        elseif strcmp(sampling,'adaptive')
            [subset_f,probs_f,var_f] = adaptive_selection_twopt(X,Xsp,Fy,lipY,batchsize,delta,old_deltas,xkin,center_ind,'trial');
        elseif strcmp(sampling,'uniform')
            [subset_f,probs_f] = uniform_selection(m,batchsize,xkin,center_ind,'trial');
        end
       
       %% Step 4.5: Update model to reflect second subset        
        F(xkin,subset_f) = fun(X(xkin,:),subset_f); 
        Eval(xkin,subset_f) = true;
        
        old_center_ind = center_ind;
        center_ind(subset_f) = xkin;

        oldFy = Fy;

        Fy(subset_f) = F(xkin,subset_f);
        Eval(xkin,subset_f) = true;
        old_deltas(subset_f) = delta;

        oldGres = Gres; oldHres = Hres; 
        
        [nf,Eval,Res,X,F,Hres,lipY,~,~,Gres] = ... 
            update_model_centers_and_models(fun,xkin,n,m,subset_f,delta,nf,Eval,center_ind,old_center_ind,Res,X,F,Fy,Gres,Hres,npmax,Par,L,U,nfmax,lipY,probs_f,printf,nolip);
       
        %% Step 4.6: Evaluate at new point
        nf = nf + 1;
        X(nf,:) = X(xkin,:)+Xsp;
        F(nf,subset_f) = fun(X(nf,:),subset_f);
        Eval(nf,subset_f) = true;

        %% Step 5: Compute estimates:
        
        [f0,fs] = ... 
            compute_estimates(oldFy,F,X,xkin,old_center_ind,oldGres,oldHres,subset_f,probs_f,subset_f,probs_f,nf);
        
        %% Step 6: Determine success
        rho = (f0-fs)/(-mdec);
        
        %% Step 7: Accept point or not
        if rho >= eta1 && (f0 - fs) > 0
            % success
            if ~first_success_now
                first_success_now = true;
            end
            xkin = nf;
            f = fs;
        else
            % not success
            f = f0;
        end
        
        %% Step 8: Trust region adjustments
        if rho >= eta1 && norm(Xsp) > 0.75*delta
            delta = min(delta*gam1,maxdelta);
        elseif rho < eta1
            delta = max(delta*gam0,mindelta);
        end               
    else % not acceptable
        delta = max(delta*gam0,mindelta);
    end % end if acceptable
    
    %% Step 8.5: Print and iterate
    if printf 
        if strcmp(sampling,'arbitrary') || strcmp(sampling,'adaptive')
            fprintf(progstr, sum(sum(Eval)), delta, f, ng, var_m, var_f); 
        else
            fprintf(progstr, sum(sum(Eval)), delta, f, ng);
        end
    end

end % end while outer loop 

end

function [subset,probs,var] = adaptive_selection_tr(X,fY,lipY,batchsize,delta,old_deltas,xkin,center_ind)
    
    pi_param = 0.01; C_param = 10*sum(lipY);
    Y = X(center_ind,:);
    x = X(xkin,:);  
    [m,n] = size(Y);
    Par2 = max(10,sqrt(n)); %hard coded, but really just a copy of Par(2) in main method. 
    
    errors = zeros(1,m);
    
    if batchsize == m
        probs = ones(1,m);
    else
        % attempt to upper-bound the errors of each model
        for comp = 1:m
            yj = Y(comp,:);
            dist = norm(yj-x);
            errors(comp) = (abs(fY(comp))*lipY(comp)/2)*...
                (3*(dist + delta)^2 + ...
                 +  sqrt(n)*Par2*(old_deltas(comp)^2*(dist+delta)) ...
                 + 3*delta^2 + sqrt(n)*Par2*delta^3);
        end

        [sorted,sortinds] = sort(errors);
        cumsorted = cumsum(sorted);
    end
       
    stopped = false; cb = batchsize;
    while ~stopped
        if m <= cb
            subset = 1:m;
            probs = ones(1,m);
            var = 0;
            break
        else
            k = m;
            while k >= 1
                if cb + k - m <= cumsorted(k)/sorted(k)
                    break
                end
                k = k - 1;
            end
            % should be unnecessary, but just in case:
            if k == 0
                k = 1;
            end
            % get optimal probabilities
            probs = ones(1,m);
            probs(sortinds(1:k)) = (cb + k - m)*errors(sortinds(1:k))/cumsorted(k);
            probs = max(probs,eps); 
            subset = get_subset(probs,cb);
        end

        var = sum(((1.0./probs)-1.0).*errors.^2);
        if var < pi_param*C_param^2*delta^4
            break
        else
            cb = cb + batchsize;
        end
    end
end

function [subset,probs,var] = adaptive_selection_twopt(X,s,fY,lipY,batchsize,delta,old_deltas,xkin,center_ind,sense)
    
    pi_param = 0.01; C_param = 10*sum(lipY);
    Y = X(center_ind,:);
    x = X(xkin,:);  
    [m,n] = size(Y);
     Par2 = max(10,sqrt(n)); %hard coded, but really just a copy of Par(2) in main method. 
    
    stepsize = norm(x-s);
    
    errors = zeros(1,m);
    
    if batchsize == m
        probs = ones(1,m);
    else
        % attempt to upper-bound the errors of each model
        for comp = 1:m
            yj = Y(comp,:);
            dist = norm(yj-x);
            errors(comp) = (abs(fY(comp))*lipY(comp)/2)*max(3*dist^2 + sqrt(n)*Par2*old_deltas(comp)^2*dist,...
                3*norm(x+s-yj)^2 + sqrt(n)*Par2*old_deltas(comp)^2*norm(x+s-yj) + 3*stepsize^2 + sqrt(n)*Par2*delta^2*stepsize);
        end

        [sorted,sortinds] = sort(errors);
        cumsorted = cumsum(sorted);
    end

    if strcmp(sense,'center')
        selectable = find(center_ind ~= xkin); ls = length(selectable);
    elseif strcmp(sense,'trial')
        selectable = 1:m; ls = m; 
    end
    
    cb = batchsize; stopped = false;
    while ~stopped
        if ls <= cb
            subset = selectable; 
            probs = ones(1,m); 
            var = 0; 
            break
        else
            k = m;
            while k >= 1
                if cb + k - m <= cumsorted(k)/sorted(k)
                    break
                end
                k = k - 1;
            end
            % should be unnecessary, but just in case:
            if k == 0
                k = 1;
            end
            % get optimal probabilities
            probs = ones(1,m);
            probs(sortinds(1:k)) = (cb + k - m)*errors(sortinds(1:k))/cumsorted(k);
            probs = max(probs,eps); 
            subset = get_subset(probs(selectable),cb);
            subset = selectable(subset);
        end

        var = sum(((1.0./probs)-1.0).*errors.^2);
        if var < pi_param*C_param^2*delta^4
            break
        else
            cb = cb + batchsize;
        end
    end
end

function [subset,probs] = uniform_selection(m,batchsize,xkin,center_ind,sense)

    probs = (batchsize/m)*ones(1,m);
    if strcmp(sense,'center')
        selectable = find(center_ind ~= xkin); ls = length(selectable);
    elseif strcmp(sense,'trial')
        selectable = 1:m; ls = m; 
    end
    if ls < batchsize
        remaining = setdiff(1:m,selectable);
        subset = union(selectable, datasample(remaining,batchsize-ls,'Replace',false,'Weights',probs(remaining)));
    else
        subset = datasample(selectable, batchsize,'Replace',false,'Weights',probs(selectable));
    end
end
