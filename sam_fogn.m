function [X,F,flag,xkin,Eval] = ...
    sam_fogn(fun,X0,n,npmax,nfmax,gtol,delta,nfs,m,F0,xkin,L,U,printf,batchsize,sampling,lipY)

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
    disp('  nf   delta           f0           g0       ');
    progstr = '%4i %9.2e  %11.5e %12.4e \n'; % Line-by-line
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

X = [X0; zeros(nfmax-1,n)]; % Stores the point locations
F = nan(nfmax,m); % Stores the function values
Eval = logical(zeros(nfmax,m)); % Stores whether eval i includes component j
nf = 1;
Fy = F(xkin,:);

Gres = zeros(n,m);
center_ind = ones(1,m); % Stores model centers

if nargin == 16
    lipY = ones(1,m); % default lipschitz constant estimates
end

% evaluate all components immediately:
crit_check = true;

while sum(sum(Eval))<nfmax && delta > mindelta
    %% Step 1: Choose a subset
    if crit_check
        subset_m = 1:m; probs_m = ones(1,m);
        subset_f = 1:m; probs_f = ones(1,m);
    else
        if strcmp(sampling,'arbitrary')
            [subset_m,probs_m] = arbitrary_selection_tr(X(center_ind,:),Fy,Gres,lipY,X(xkin,:),batchsize,delta);
        elseif strcmp(sampling,'uniform')
            [subset_m,probs_m] = uniform_selection(m,batchsize);
        else
            error('sampling must be arbitrary or uniform');
        end
    end
    
    %% Step 2: Update model centers and models
    oldGres = Gres;
    
    old_center_ind = center_ind;
    center_ind(subset_m) = xkin;
    [F(xkin,subset_m),Gres(:,subset_m)] = fun(X(xkin,:),subset_m); % PRACTICAL: THIS MAY BE LIKELY EVALUATING REDUNDANT COMPONENTS
    
    oldFy = Fy;    
    Fy(subset_m) = F(xkin,subset_m); 
    
    Eval(xkin,subset_m) = true;
    
    G = zeros(n,1); 
    H = zeros(n); 
    oldGres = Gres;
    for j = 1:m
       
        Dj = X(xkin,:)-X(old_center_ind(j),:); 
        old_cj = oldFy(j) + Dj*oldGres(:,j); %+ 0.5*Dj*Hres(:,:,j)*Dj';
        old_gj = oldGres(:,j); %+ Hres(:,:,j)*Dj';

        G = G + old_cj*old_gj; 
        H = H + 2*oldGres(:,j)*oldGres(:,j)'; %+ 2*old_cj*Hres(:,:,j); 
        if ismember(j,subset_m)
            % add in the new model information
            new_cj = Fy(j);
            new_gj = Gres(:,j); 
            G = G - (1.0/probs_m(j))*old_cj*old_gj + (1.0/probs_m(j))*new_cj*new_gj;
            Hj = - (2.0/probs_m(j))*((oldGres(:,j)*oldGres(:,j)')) ...
                + (2.0/probs_m(j))*((Gres(:,j)*Gres(:,j)'));  
            H = H + Hj;
            
            % update Lipschitz estimates
            newlip = norm(oldGres(:,j)-Gres(:,j))/norm(Dj);
            if norm(Dj) > 0 && newlip > lipY(j)
                lipY(j) = newlip;
            end
        end       
    end
        
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
            [subset_f,probs_f] = arbitrary_selection_twopt(X(center_ind,:),Fy,Gres,lipY,X(xkin,:),batchsize,Xsp);
        elseif strcmp(sampling,'uniform')
            [subset_f,probs_f] = uniform_selection(m,batchsize);
        end
       
       %% Step 4.5: Update model to reflect second subset
        center_ind(subset_f) = xkin;
        
        % Do some function evaluations
        % Always compute at the trial step for now:
        nf = nf + 1;
        X(nf,:) = X(xkin,:)+Xsp;
        F(nf,subset_f) = fun(X(nf,:),subset_f);
        Eval(nf,subset_f) = true;
        % And then evaluate the center point if necessary:
        for j = subset_f % do "in parallel"
           if ~Eval(xkin,j)
              F(xkin,j) = fun(X(xkin,:),j); 
              Eval(xkin,j) = true;
           end
        end

        %% Step 5: Compute estimates:
        f0 = 0; fs = 0;
        for j = 1:m
            % f0: 
            if isnan(oldFy(j)) % this was the first time we ever evaluated it. 
                f0 = f0 + F(xkin,j)^2;
            else
                Dj = X(xkin,:)-X(old_center_ind(j),:); 
                old_cj = oldFy(j) + Dj*oldGres(:,j); %+ 0.5*Dj*oldHres(:,:,j)*Dj';
                f0 = f0 + old_cj^2;
                if ismember(j,subset_f)
                    f0 = f0 + (1.0/probs_f(j))*(F(xkin,j)^2 - old_cj^2);
                end
            end 
            % fs: 
            if isnan(oldFy(j))
                fs = fs + F(nf,j)^2;
            else
                Dj = X(nf,:)-X(old_center_ind(j),:); 
                old_cj = oldFy(j) + Dj*oldGres(:,j); %+ 0.5*Dj*oldHres(:,:,j)*Dj'; %fogn
                fs = fs + old_cj^2;
                if ismember(j,subset_f)
                    fs = fs + (1.0/probs_f(j))*(F(nf,j)^2 - old_cj^2);
                end
            end
        end
        
        %% Step 6: Determine success
        rho = (f0-fs)/(-mdec);
        
        %% Step 7: Accept point or not
        if rho >= eta1 && (f0 - fs) > 0
            % success
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
    if printf; fprintf(progstr, sum(sum(Eval)), delta, f, ng); end

end % end while outer loop 

end

function [subset,probs] = arbitrary_selection_tr(Y,fY,Gres,lipY,x,batchsize,delta)
    
    m = size(Y,1);
    
    errors = zeros(1,m);
    
    if batchsize == m
        probs = ones(1,m);
    else
        % attempt to upper-bound the errors of each model
        for comp = 1:m
            yj = Y(comp,:);
            dist = norm(yj-x);
            Mj = abs(fY(comp)) + abs((yj-x)*Gres(:,comp)) + (lipY(comp)/2)*dist^2;
            errors(comp) = (abs(fY(comp))*lipY(comp)/2)*(delta^2 + (dist+delta)^2) ...
                + (lipY(comp)*Mj)/2*delta^2;
        end

        [sorted,sortinds] = sort(errors);
        cumsorted = cumsum(sorted);
        k = m;
        while k >= 1
            if batchsize + k - m <= cumsorted(k)/sorted(k)
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
        probs(sortinds(1:k)) = (batchsize + k - m)*errors(sortinds(1:k))/cumsorted(k);
    end
    subset = datasample(1:m,batchsize,'Replace',false,'Weights',probs);
end

function [subset,probs] = arbitrary_selection_twopt(Y,fY,Gres,lipY,x,batchsize,s)
    
    m = size(Y,1);
    
    errors = zeros(1,m);
    
    stepsize = norm(x-s);
    
    if batchsize == m
        probs = ones(1,m);
    else
        % attempt to upper-bound the errors of each model
        for comp = 1:m
            yj = Y(comp,:);
            dist = norm(yj-x);
            Mj = abs(fY(comp)) + abs((yj-x)*Gres(:,comp)) + (lipY(comp)/2)*dist^2;
            errors(comp) = max((abs(fY(comp))*lipY(comp)/2)*dist^2,...
                (abs(fY(comp))*lipY(comp)/2)*norm(x+s-yj)^2 + ...
                (lipY(comp)*Mj/2)*stepsize^2);
        end

        [sorted,sortinds] = sort(errors);
        cumsorted = cumsum(sorted);
        k = m;
        while k >= 1
            if batchsize + k - m <= cumsorted(k)/sorted(k)
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
        probs(sortinds(1:k)) = (batchsize + k - m)*errors(sortinds(1:k))/cumsorted(k);
    end
    try
        subset = datasample(1:m,batchsize,'Replace',false,'Weights',probs);
    catch
        % in case the errors are zero at all points 
        probs = (batchsize/m)*ones(1,m);
        subset = datasample(1:m,batchsize,'Replace',false,'Weights',probs);
    end
end

function [subset,probs] = uniform_selection(m,batchsize)

    probs = (batchsize/m)*ones(1,m);
    subset = datasample(1:m,batchsize,'Replace',false,'Weights',probs);

end

% function [subset,probs] = adaptive_selection(Y,lipY,x,delta)
%     
%     [m,n] = size(Y);
%     
%     errors = zeros(1,m);
%     
%     batchsize = 1; keepgoing = true;
%     alpha = 0.95; C = m*max(lipY);
%     
%     % upper-bound the errors of each model
%     for comp = 1:m
%         yj = Y(comp,:);
%         dist = norm(yj-x);
%         errors(comp) = (0.5*lipY(comp))*(dist^2 + 2*delta*dist + 2*delta^2);
%     end
%     [sorted,sortinds] = sort(errors);
%     cumsorted = cumsum(sorted);
%     
%     while keepgoing
%         if batchsize == m
%             probs = ones(1,m);
%         else           
%             k = m;
%             while k >= 1
%                 if batchsize + k - m <= cumsorted(k)/sorted(k)
%                     break
%                 end
%                 k = k - 1;
%             end
%             % should be unnecessary, but just in case:
%             if k == 0
%                 k = 1;
%             end
%             % get optimal probabilities
%             probs = ones(1,m);
%             probs(sortinds(1:k)) = (batchsize + k - m)*errors(sortinds(1:k))/cumsorted(k);
%         end
%         
%         % ready to stop yet?
%         v = sum(((1.0./probs) - 1.0).*errors);
%         if v <= (1-alpha)*C^2*delta^4
%             keepgoing = false;
%         else
%             batchsize = batchsize + 1;
%         end
%     end
%     try
%         subset = datasample(1:m,batchsize,'Replace',false,'Weights',probs);
%     catch
%         [subset,probs] = uniform_selection(m,batchsize);
%     end
% end