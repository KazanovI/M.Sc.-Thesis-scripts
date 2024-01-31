function [cells_peak] = cell_activity_stat(dat,border)
%cell_activity_stat get window(border) activity and check if statistically
%significant
    base = cellfun(@(x) squeeze(mean(x(:,border(1,1):border(1,2),:),2,'omitnan')),dat,'UniformOutput',false);
    resp = cellfun(@(x) squeeze(mean(x(:,border(2,1):border(2,2),:),2,'omitnan')),dat,'UniformOutput',false);
    cells_peak = cell(1,length(base));
    for j = 1:length(base) % session
        cells_peak{j} = zeros(size(base{j},1),3);
        for c = 1:size(base{j},1) % per cell
            ind_p = find(~isnan(base{j}(c,:)));
            if ~isempty(ind_p)
                %  Wilcoxon signed rank test
                comp_c = [base{j}(c,ind_p)' resp{j}(c,ind_p)'];
                % Enhanced
                [~,h1] = signrank(comp_c(:,1),comp_c(:,2),'tail','left','alpha',0.025);
                cells_peak{j}(c,1) = h1;
                % Suppressed
                [~,h2] = signrank(comp_c(:,1),comp_c(:,2),'tail','right','alpha',0.025);
                cells_peak{j}(c,2) = h2;
                % Neutral (No change)
                if ~any(cells_peak{j}(c,1:2))
                    cells_peak{j}(c,3) = 1;
                end
            end
        end           
    end
end