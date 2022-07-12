function [hlu,hld] = plot_trajectory_percentiles2(H,prob_min,uniform_solvers,dynamic_solvers,width,mode,batchsizes,m)

% H - nf x np x ns array generated by run
% prob_min - np length array of minimum problem data
% solvers_to_plot - indices (3rd dim of H) to show in plot
colors=get(gca,'colororder');
markers={'s','o','^','v','p','<','x','h','+','d','*'};
LW=2;
FS=16;
MS=8;
transparency = 0.2; 

% for clear plots:
%H((H<2e-32)) = 2e-32;
[nf,~,~] = size(H); % Grab the dimensions

%% Process/plot the uniform solvers first
% Produce a suitable history array with sorted entries:
for j = uniform_solvers
    for i = 2:nf
        H(i,:,j) = min(H(i,:,j),H(i-1,:,j));
    end
end

l = 1;
ns = length(uniform_solvers);
for s = uniform_solvers % each solver     
    data = H(:,:,s)'-repmat(prob_min',1,nf);
    low = prctile(data,50-width);
    med = prctile(data,50);
    upp = prctile(data,50+width);
    
    for i = 2:nf
        if med(i) == 0
            med(i) = med(i-1);
        end
    end
    
    if strcmp(mode,'eval')
        xaxis = (1:nf)/m;
    elseif strcmp(mode,'batch')
        xaxis = (1:nf)/batchsizes(l);
    end
    hlu(l) = semilogy(xaxis,med,'LineWidth',LW,'LineStyle','-','Color',colors(mod(l-1,ns)+1,:),...
        'MarkerSize',MS,'Marker',markers{mod(l-1,ns)+1},'MarkerIndices',floor(linspace(1,nf,10)),...
        'MarkerFaceColor',colors(mod(l-1,ns)+1,:),'MarkerEdgeColor',colors(mod(l-1,ns)+1,:));
    hold on
    
    if width > 0
        for i = 2:nf
            if low(i) == 0
                low(i) = low(i-1);
            end
            if upp(i) == 0
                upp(i) = upp(i-1);
            end
        end
        hlcu(l) = fill([xaxis,fliplr(xaxis)],[low,fliplr(upp)],colors(mod(l-1,ns)+1,:));
        set(hlcu(l),'FaceAlpha',transparency,'EdgeColor',colors(mod(l-1,ns)+1,:));
    end

    % lines on the boundaries - no markers
    semilogy(xaxis,low,'LineWidth',LW-1,'LineStyle','-','Color',colors(mod(l-1,ns)+1,:));
    semilogy(xaxis,upp,'LineWidth',LW-1,'LineStyle','-','Color',colors(mod(l-1,ns)+1,:));
    l = l + 1;

end

%% Now do dynamic solvers:

for j = dynamic_solvers
    for i = 2:nf
         H(i,:,j) = min(H(i,:,j),H(i-1,:,j));
    end
end

l = 1;
for s = dynamic_solvers % each solver     
    data = H(:,:,s)'-repmat(prob_min',1,nf);
    low = prctile(data,25);
    med = prctile(data,50);
    upp = prctile(data,75);
    
    for i = 2:nf
        if med(i) == 0
            med(i) = med(i-1);
        end
    end
    if strcmp(mode,'eval')
        xaxis = (1:nf)/m;
    elseif strcmp(mode,'batch')
        xaxis = (1:nf)/batchsizes(l);
    end
    hld(l) = semilogy(xaxis,med,'LineWidth',LW,'LineStyle','--','Color',colors(mod(l-1,ns)+1,:),...
        'MarkerSize',MS,'Marker',markers{mod(l-1,ns)+1},'MarkerIndices',floor(linspace(1,nf,10)),...
        'MarkerFaceColor',colors(mod(l-1,ns)+1,:),'MarkerEdgeColor',colors(mod(l-1,ns)+1,:));
    hold on
    
    if width > 0
        for i = 2:nf
            if isnan(low(i)) || low(i) == 0
                low(i) = low(i-1);
            end
            if isnan(upp(i)) || upp(i) == 0
                upp(i) = upp(i-1);
            end
        end
        hlcd(l) = fill([xaxis,fliplr(xaxis)],[low,fliplr(upp)],colors(mod(l-1,ns)+1,:));
        set(hlcd(l),'FaceAlpha',transparency,'EdgeColor',colors(mod(l-1,ns)+1,:));
    end
    
    semilogy(xaxis,low,'LineWidth',LW-1,'LineStyle','--','Color',colors(mod(l-1,ns)+1,:));
    semilogy(xaxis,upp,'LineWidth',LW-1,'LineStyle','--','Color',colors(mod(l-1,ns)+1,:));
    l = l + 1;
end

hold on

ax = gca;
xlim(ax,[0,nf/m]);
%ylim(ax,[yl_min,yl_max])
%set(ax,'YTick',linspace(yl_min,yl_max,10))




end