function [f0,fs] = compute_estimates(oldFy,F,X,xkin,old_center_ind,oldGres,oldHres,subset_fc,probs_fc,subset_fs,probs_fs,nf)

m = length(oldFy);

f0 = 0; fs = 0;

for j = 1:m
    % f0: 
    if isnan(oldFy(j)) % this was the first time we ever evaluated it. 
        f0 = f0 + F(xkin,j)^2;
    else
        Dj = X(xkin,:)-X(old_center_ind(j),:); 
        old_cj = oldFy(j) + Dj*oldGres(:,j) + 0.5*Dj*oldHres(:,:,j)*Dj';
        f0 = f0 + old_cj^2;
        if ismember(j,subset_fc)
            f0 = f0 + (1.0/probs_fc(j))*(F(xkin,j)^2 - old_cj^2);
        end
    end 
    % fs: 
    if isnan(oldFy(j))
        fs = fs + F(nf,j)^2;
    else
        Dj = X(nf,:)-X(old_center_ind(j),:); 
        old_cj = oldFy(j) + Dj*oldGres(:,j) + 0.5*Dj*oldHres(:,:,j)*Dj';
        fs = fs + old_cj^2;
        if ismember(j,subset_fs)
            fs = fs + (1.0/probs_fs(j))*(F(nf,j)^2 - old_cj^2);
        end
    end
end

