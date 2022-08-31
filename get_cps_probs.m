function [pi_tilde,lambda,jip] = get_cps_probs(pi,b)

% COMPUTES DISTRIBUTION FOR CONDITIONAL POISSON SAMPLING THAT MATCHES THE
% INCLUSION PROBABILITIES FOR POISSON SAMPLING 

% initial point
pi_tilde = pi; p = length(pi_tilde);

residual = inf; maxiters = 1000; iter = 0;
while residual > p*sqrt(eps) && iter < maxiters
    psi = recursive(pi_tilde,b);
    %Delta = pi - pi.^2;
    %grad = Delta./(pi_tilde.*(1-pi_tilde));
    alpha = 1.0; residual_old = norm(pi-psi); residual = inf;
    while residual >= residual_old && alpha > eps
        pi_tilde_new = pi_tilde + alpha*(pi-psi);
        %pi_tilde_new = pi_tilde + (pi - psi);
        residual = norm(pi-psi);
        alpha = 0.5*alpha;
    end
    pi_tilde = pi_tilde_new;
    iter = iter + 1;
end

% run nonlinear least squares
% fun = @(pi_tilde)fitting_function(pi_tilde,b,pi);
% %options = optimoptions('lsqnonlin','Display','off','OptimalityTolerance',1e-16);
% %[pi_tilde,residual_norm] = lsqnonlin(fun,pi_tilde,zeros(p,1),ones(p,1),options);
% options = optimoptions('fmincon','Display','off','OptimalityTolerance',1e-16);
% [pi_tilde] = fmincon(fun,pi_tilde,[],[],ones(1,p),b,zeros(p,1),ones(p,1),[],options);

% now pick a scalar c so that pi_tilde sums to b: 
lambda = log(pi_tilde./(1.0-pi_tilde));

% golden section search
% lb = -b;
% ub = b;
% gr = (sqrt(5) + 1)/2;
% m1 = ub - (ub-lb)/gr;
% m2 = lb + (ub-lb)/gr;
fun = @(c)totalnewpi(c,lambda,b);
% while abs(ub-lb) > sqrt(eps)
%     if fun(m1) < fun(m2)
%         ub = m2;
%     else
%         lb = m1;
%     end
% 
%     m1 = ub - (ub-lb)/gr;
%     m2 = lb + (ub-lb)/gr;
% end
% c = (ub+lb)/2;

options = optimoptions('fminunc','Display','off');
try
    c = fminunc(fun,0,options);
catch
    pause()
end
[~,pi_tilde] = totalnewpi(c,lambda,b);
%residual_norm = sum(fitting_function(pi_tilde,b,pi).^2);

psi = recursive(pi_tilde,b); lambda = lambda + c;

if nargout == 3 % compute joint inclusion probs
jip = zeros(p); dolater_k = []; dolater_l = []; 
for k = 1:p
    for l = (k+1):p
        if abs(lambda(k) - lambda(l)) > sqrt(eps)
            jip(k,l) = (psi(k)*exp(lambda(l)) - psi(l)*exp(lambda(k)))/(exp(lambda(l))-exp(lambda(k)));
        else
            dolater_k = cat(1,dolater_k,k); dolater_l = cat(1,dolater_l,l);
        end
    end
end
jip = jip + jip'; jip = jip + diag(psi);
remaining = length(dolater_k);
if remaining > 0
    for r = 1:remaining
        k = dolater_k(r); l = dolater_l(r);
        inds = find(abs(lambda - lambda(l)) <= sqrt(eps));
        skippedinds = setdiff(1:p,inds);
        sumnoneq = sum(jip(k,skippedinds));
        jip(k,l) = 1/(length(inds)-1)*(psi(l)*(b-1) - sumnoneq);
        jip(l,k) = jip(k,l);
    end
end
end
end

function [psi] = recursive(pi_tilde,b)

lower_psi = zeros(1,length(pi_tilde));
for k = 1:b    
    psi = (pi_tilde./(1.0-pi_tilde)).*(1.0 - lower_psi);

    lower_psi = k*psi/sum(psi);
end

psi = lower_psi;

end

% function fval = fitting_function(pi_tilde,b,pi)
% 
% residual = recursive(pi_tilde,b) - pi';
% 
% fval = sum(residual.^2);
% 
% end

function [f,newpi] = totalnewpi(c,lambda,b)
    newpi = exp(lambda + c)./(1.0 + exp(lambda+c));
    f = (sum(newpi)-b)^2;
end