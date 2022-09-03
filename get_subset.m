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
        new_probs = get_cps_probs(smaller_probs,b);
        if abs(sum(new_probs) - b) > 1.0
            % THIS IS A SAFEGUARD IF SOMETHING WENT NUMERICALLY AWRY. IT
            % SHOULD TRIGGER RARELY.  
            % 1.0 IS ARBITRARY, BUT INDICATES SOMETHING HAS GONE QUITE
            % WRONG. 
            warning('Something went far too wrong in the CPS computation. Doing a greedy thing on this iteration.')
            [~,new_subset_inds] = maxk(new_probs,b); 
            new_subset = notinds(new_subset_inds);
            rejected = false;
        else
            rejected = true;
        end
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
end