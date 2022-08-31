function subset = get_subset(probs,b)

    p = length(probs);

    % FIRST, GET PROBS FOR CONDITIONAL POISSON SAMPLING SUCH THAT
    % FIRST-ORDER INCLUSION PROBABILITIES MATCH

    % GET RID OF PROB ONES:
    inds = find(probs>=1.0-eps);
    subset = inds; 
    notinds = setdiff(1:p,inds);
    smaller_probs = probs(notinds);
    b = b - length(subset);
    if b > 0
         [new_probs,lambda] = get_cps_probs(smaller_probs,b);
%         if abs(sum(new_probs) - b) > 0.0001
%             pause()
%         end
%         newinds = [];
%         for k = 1:b
%             if k == 1
%                 q = new_probs/b;
%             elseif k >= 2
%                 % update q
%                 qnewind = q(newind);
%                 equal = find(q==qnewind);
%                 unequal = setdiff(1:length(notinds),equal);
%                 
%                 sumq = 0;
%                 for j = 1:length(notinds)
%                     if ismember(j,unequal) && ~ismember(j,newinds)
%                         elnewind = exp(lambda(newind));
%                         eldiff = exp(lambda(j)-lambda(newind));
%                         num = elnewind*((q(j)/qnewind)-eldiff);
%                         denom = elnewind*(1.0-eldiff);
%                         tmp = (1/(b-(k-2)-1))*(num/denom);
%                         q(j) = tmp;
%                         sumq = sumq + q(j);
%                     elseif ismember(j,newinds)
%                         q(j) = 0;
%                     end
%                 end
%                 card = sum(ismember(1:length(notinds),equal).*(~ismember(1:length(notinds),newinds)));
%                 for j = 1:length(notinds)
%                     if ismember(j,equal) && ~ismember(j,newinds)
%                         q(j) = (1 - sumq)/card;
%                     end
%                 end
%             end
% %             if abs(sum(q)-1.0) > sqrt(eps)
% %                 pause();
% %             end
%             % repairing small errors: 
%             q(q<0) = 0;
%             q = q/sum(q);          
%             try
%                 newind = datasample(1:length(notinds),1,'Weights',q);
%             catch
%                 % This should not trigger very often - sometimes, when some
%                 % q become very close to 0 (typically near the end of an
%                 % optimization run), numerical problems happen and q
%                 % becomes NaN - this catch addresses that by doing a less
%                 % disciplined random draw
%                 %[~,new_subset] = maxk(new_probs,b);
%                 new_subset = datasample(1:length(notinds),b,'Weights',smaller_probs);
%                 subset = [inds notinds(new_subset)];
%                 subset = unique(subset);
%                 %warning('Catch happened in get_subset')
%                 break
%             end
%             if isempty(subset)
%                 subset = notinds(newind);
%                 newinds = newind;
%             else
%                 newinds = cat(2,newinds,newind);
%                 subset = cat(2,subset,notinds(newind));
%             end
%         end
%     end  
        rejected = true;
        while rejected
            % POISSON SAMPLE
            new_subset = [];
            for k = 1:length(notinds)
                % flip coin
                if rand() < new_probs(k)
                    new_subset = cat(1,new_subset,notinds(k));
                end
            end
            if length(new_subset) == b
                rejected = false;
            end
        end
        if isempty(subset)
            subset = new_subset(:);
        else
            subset = cat(1,subset(:),new_subset(:));
        end
    end
% end