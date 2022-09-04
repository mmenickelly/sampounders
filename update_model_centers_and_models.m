function [nf,Eval,Res,X,F,Hres,lipY,G,H,Gres,Lip] = update_model_centers_and_models(fun,xkin,n,m,subset,delta,nf,Eval,center_ind,old_center_ind,Res,X,F,Fy,Gres,Hres,npmax,Par,L,U,nfmax,lipY,probs,printf,nolip,Lip) 
       
Hresdel = zeros(n,n,m); %Gres=zeros(n,m);
valid = zeros(1,m); np = zeros(1,m);
for j = subset(:)'
    yj = old_center_ind(j);
    for i=1:nf
        if Eval(i,j)
            D = X(i,:)-X(yj,:);
            Res(i,j) = (F(i,j)-Fy(j))-.5*D*Hres(:,:,j)*D';
        end
    end
    %Check all models, for now make sure we have np=n for all
    %! Warning, what follows assumes one model at a time, ideally want
    % batches.
    [Mdira,np(j)] = ...
        formquad_indep(X(Eval(:,j),:),Res(Eval(:,j),j),delta,X(yj,:),npmax,Par,0);
    if np(j)<n  % Must obtain and evaluate bounded geometry points
        [Mdira,np(j)] = bmpts(X(yj,:),Mdira(1:n-np(j),:),L,U,delta,Par(3));
        for i=1:min(n-np(j),nfmax-nf)
            nf = nf+1;
            X(nf,:) = min(U,max(L,X(yj,:)+Mdira(i,:))); % Temp safeguard
            F(nf,j) = fun(X(nf,:),j);
            Eval(nf,j) = true;
            if printf
                fprintf('Geometry point evaluated for model %i\n',j);
            end
            D = Mdira(i,:);
            Res(nf,j) = (F(nf,j)-Fy(j))-.5*D*Hres(:,:,j)*D';
        end
        if nf>=nfmax; break; end
    end
    [~,np(j),valid(j),Gres(:,j),Hresdel(:,:,j)] = ...
        formquad_indep(X(Eval(:,j),:),Res(Eval(:,j),j),delta,X(yj,:),npmax,Par,0);
end

Hres = Hres + Hresdel;
G = zeros(n,1); H = zeros(n); Hresdel = zeros(n,n,m);
oldGres = Gres; oldHres = Hres; 
for j = 1:m
    Dj = X(xkin,:)-X(old_center_ind(j),:); 
    old_cj = Fy(j) + Dj*Gres(:,j) + 0.5*Dj*Hres(:,:,j)*Dj';
    old_gj = Gres(:,j) + Hres(:,:,j)*Dj';
    G = G + old_cj*old_gj;
    H = H + 2*old_cj*Hres(:,:,j) + 2*Gres(:,j)*Gres(:,j)';
    if ismember(j,subset)
        % also need to compute model at the NEW index
        oldGresj = Gres(:,j);
        yj = center_ind(j);
        for i=1:nf
            if Eval(i,j)
                D = X(i,:)-X(yj,:);            
                Res(i,j) = (F(i,j)-Fy(j))-.5*D*Hres(:,:,j)*D';
            end
        end
        %Check all models, for now make sure we have np=n for all
        %! Warning, what follows assumes one model at a time, ideally want
        % batches.
        [Mdira,npa] = ...
            formquad_indep(X(Eval(:,j),:),Res(Eval(:,j),j),delta,X(yj,:),npmax,Par,0);
        if npa<n  % Must obtain and evaluate bounded geometry points
            [Mdira,npa] = bmpts(X(yj,:),Mdira(1:n-npa,:),L,U,delta,Par(3));
            for i=1:min(n-npa,nfmax-nf)
                nf = nf+1;
                X(nf,:) = min(U,max(L,X(yj,:)+Mdira(i,:))); % Temp safeguard
                F(nf,j) = fun(X(nf,:),j);
                Eval(nf,j) = true;
                if printf
                    fprintf('Geometry point evaluated for model %i\n',j);
                end
                D = Mdira(i,:);
                Res(nf,j) = (F(nf,j)-Fy(j))-.5*D*Hres(:,:,j)*D';
            end
            if nf>=nfmax; break; end
        end
        [~,np(j),valid(j),Gres(:,j),Hresdel(:,:,j)] = ...
                formquad_indep(X(Eval(:,j),:),Res(Eval(:,j),j),delta,X(yj,:),npmax,Par,0);
        % add in the new model information
        new_cj = Fy(j);
        new_gj = Gres(:,j); 
        G = G - (1.0/probs(j))*old_cj*old_gj + (1.0/probs(j))*new_cj*new_gj;
        Hj = - (2.0/probs(j))*(old_cj*Hres(:,:,j) + (oldGresj*oldGresj')) ...
            + (2.0/probs(j))*(new_cj*(Hres(:,:,j)+Hresdel(:,:,j)) + (Gres(:,j)*Gres(:,j)'));  
        H = H + Hj;

        % update lipschitz estimates
        newlip = norm(oldGresj-Gres(:,j))/norm(Dj);
        if norm(Dj) > 0 %&& nolip && newlip > lipY(j)
            lipY(j) = newlip;
            Lip(xkin,j) = newlip;
        end
    end       
end
Hres = Hres + Hresdel;