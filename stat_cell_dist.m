function stat_cell_dist(stat_inds,locs)
    clrs = get_color(["Enhanced","Suppressed","None"]);
    dists = cell(1,length(stat_inds));
    for s_i = 1:length(stat_inds)
        dists{s_i} = cell(1,size(stat_inds{s_i},2));
        for ct = 1:size(stat_inds{s_i},2)
            ct_locs = locs{s_i}(logical(stat_inds{s_i}(:,ct)),:);
            inds = cell(size(ct_locs,1));
            for i = 1:size(inds,1)
                for j = 1:size(inds,1)
                    inds{i,j} = [i j];
                end
            end
            inds_ind = tril(ones(size(inds,1)),-1);
            inds_mat = cell2mat(inds(inds_ind ~= 0));
            dists{s_i}{ct} = zeros(size(inds_mat,1),1);
            for cell_i = 1:size(inds_mat,1)
                dists{s_i}{ct}(cell_i) = sqrt(diff(ct_locs(inds_mat(cell_i,:),1))^2 ...
                + diff(ct_locs(inds_mat(cell_i,:),2))^2);
            end 
        end
%         figure();
%         for ct = 1:size(stat_inds{s_i},2)
%             plot(locs{s_i}(logical(stat_inds{s_i}(:,ct)),1),...
%                 locs{s_i}(logical(stat_inds{s_i}(:,ct)),2),'o','Color',clrs(ct,:),'MarkerFaceColor',clrs(ct,:));
%             hold on
%         end
%         box off
%         close;
    end
    [m_dists,std_dists] = deal(cell(1,length(dists)));
    for s_i = 1:length(dists)
        m_dists{s_i} = cellfun(@(x) mean(x),dists{s_i});
        std_dists{s_i} = cellfun(@(x) std(x)/sqrt(length(x)),dists{s_i});
    end
    m_all = mean(cell2mat(m_dists'),'omitnan');
    std_all = mean(cell2mat(std_dists'),'omitnan');
    figure();
    bb = bar(m_all');
    bb.FaceColor = 'flat';
    bb.CData = clrs;
    hold on
    errorbar(m_all,std_all,'k','LineStyle','none');
    box off
    xticklabels(["Enhanced","Suppressed","None"]);
    ylabel('Distance');
    % anova1
    all_dists = cell(1,3);
    for i = 1:length(dists)
        for j = 1:length(dists{i})
            all_dists{j} = [all_dists{j} ; dists{i}{j}];
        end
    end
    num_pw = cellfun(@length,all_dists);
    lbls = repelem([1 2 3]',length(m_dists));
    dat = cell2mat(m_dists');
    dat = dat(:);
    nan_ind = isnan(dat);
    dat(nan_ind) = [];
    lbls(nan_ind) = [];
    [p,tbl,stats] = anova1(dat,lbls);
    c = multcompare(stats);

end