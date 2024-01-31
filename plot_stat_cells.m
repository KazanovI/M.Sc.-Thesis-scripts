function [perc_ret] = plot_stat_cells(cell_stats,sess_info,phase)
    global num_grp 
    if phase ~= "mix";
        colors = get_color("Enhanced","Suppressed","None");
    else
        colors = get_color("mix Enhanced","mix Suppressed","mix None");
    end
    % how many stat cells in each bin
    cells_bin = cell(size(num_grp,1),1);
    cells_bin_sess = cell(size(num_grp,1),1);
    for b = 1:size(num_grp,1) % per bin
        bin_count = zeros(1,3);
        for s_i = 1:length(cell_stats)% session
            if ismember(str2num(sess_info.Session(s_i)),num_grp(b,1):num_grp(b,2))
                bin_count = bin_count + sum(cell_stats{s_i});
                cells_bin_sess{b} = [cells_bin_sess{b}; (sum(cell_stats{s_i})/length(cell_stats{s_i}))*100];;
            end
        end
        cells_bin{b} = bin_count;
    end
    perc_ret = cells_bin_sess';
    %% anova cells perc - late addition
%     perc_ret_data = cell2mat(cellfun(@(x) x(:),perc_ret,'UniformOutput',false)');
%     bin_ind = repelem(1:6',cellfun(@length, perc_ret).*3)';
%     cell_ind = repelem(repmat([1 2 3],1,6)',repmat(cellfun(@length, perc_ret),1,3));
%     [~,~,stats] = anovan(perc_ret_data,{cell_ind bin_ind},'model','interaction','varnames',{'Cell type','Bin'});
%     [c,~,~,nms] = multcompare(stats,'Dimension',[1 2]);
%     % only enh sup
%     enh_sup = cellfun(@(x) x(:,1:2),perc_ret,'UniformOutput',false);
%     enh_sup_data = cell2mat(cellfun(@(x) x(:),enh_sup,'UniformOutput',false)');
%     bin_ind = repelem(1:6',cellfun(@length, enh_sup).*2)';
%     cell_ind = repelem(repmat([1 2],1,6)',repmat(cellfun(@length, enh_sup),1,2));
%    [~,~,stats] = anovan(enh_sup_data,{cell_ind bin_ind},'model','interaction','varnames',{'Cell type','Bin'});
%    [c,~,~,nms] = multcompare(stats,'Dimension',[1 2]);
%    % check only 1 way enh
%    enh_data = cellfun(@(x) x(:,1),perc_ret,'UniformOutput',false);
%    enh_data_vec = cell2mat(enh_data');
%    enh_bin_ind = repelem(1:6',cellfun(@length, enh_data))';
%    anova1(enh_data_vec,enh_bin_ind)
    %%
    figure();
    data = cell2mat(cells_bin);
    data_perc = (data./sum(data,2))*100;
    b=bar(data_perc);
    b(1).FaceColor = colors(1,:);
    b(2).FaceColor = colors(2,:);
    b(3).FaceColor = colors(3,:);
    legend("Enhanced","Suppressed","None");
%     title("% of cells in each cell type");
%     ylabel("%");
%     xlabel("Bin");
    box off
    
    % percentage of each cell type by bin    
    figure();
    for b = 1:size(num_grp,1)
        errorbar(b,mean(cells_bin_sess{b}(:,1)),std(cells_bin_sess{b}(:,1))/sqrt(length(cells_bin_sess{b})),...
            'Marker','+','Color',colors(1,:));
        hold on
        errorbar(b,mean(cells_bin_sess{b}(:,2)),std(cells_bin_sess{b}(:,2))/sqrt(length(cells_bin_sess{b})),...
            'Marker','+','Color',colors(2,:));
        ss1 = scatter(repmat(b-0.2,1,length(cells_bin_sess{b})),cells_bin_sess{b}(:,1),...
            'filled','MarkerFaceAlpha',0.2, 'jitter','on', 'jitterAmount',0.1);
        ss1.MarkerFaceColor = colors(1,:);
        ss2 = scatter(repmat(b+0.2,1,length(cells_bin_sess{b})),cells_bin_sess{b}(:,2),...
            'filled','MarkerFaceAlpha',0.2, 'jitter','on', 'jitterAmount',0.1);
        ss2.MarkerFaceColor = colors(2,:);
    end
    ylim([0 100]);
    xlim([0 size(num_grp,1)+1]);
    ylim([0 80]);
%     trs = find(num_grp(:,2)<=7);
%     tre = setdiff(1:size(num_grp,1),trs);
%     titl = repelem("TR",1,length(trs))';
%     titl = append(titl,string(trs));
%     titl = [titl ; append(repelem("TE",length(tre),1),string(tre'))];
    xticks(1:size(num_grp,1));
%     xticklabels(titl);
%     xlabel("Bin");
%     ylabel("percentage of cells (%)");
    box off
    % check corr
    perc_mat = cell2mat(cellfun(@(x) mean(x),cells_bin_sess,'UniformOutput',false));
    [r1,p1] = corr((1:size(perc_mat,1))',perc_mat(:,1));
    [r2,p2] = corr((1:size(perc_mat,1))',perc_mat(:,2));
%     title({["cells grouped by activity during training"];...
%         ['Corr Enh: ' num2str(r1) ', p-value: ' num2str(p1)];...
%         ['Corr Supp: ' num2str(r2) ', p-value: ' num2str(p2)]});
    enh_pnts_fit = polyfit(1:size(perc_mat,1),perc_mat(:,1),1);
    enh_line = polyval(enh_pnts_fit,1:size(perc_mat,1));
    sup_pnts_fit = polyfit(1:size(perc_mat,1),perc_mat(:,2),1);
    sup_line = polyval(sup_pnts_fit,1:size(perc_mat,1));
    plot(1:length(enh_line),enh_line,'--r',1:length(sup_line),sup_line,'--b');
    legend([ss1 ss2],"Enhanced","Suppressed");
    % inset for none cells
    figure();
    for b = 1:size(num_grp,1)
        errorbar(b,mean(cells_bin_sess{b}(:,3)),std(cells_bin_sess{b}(:,3))/sqrt(length(cells_bin_sess{b})),'Color',[0.4660 0.6740 0.1880]...
        ,'Marker','+','Color',colors(3,:),'LineStyle','none');
        hold on
        ss1 = scatter(repmat(b,1,length(cells_bin_sess{b})),cells_bin_sess{b}(:,3),[]...
            ,'filled','MarkerFaceAlpha',0.2, 'jitter','on', 'jitterAmount',0.1);
        ss1.MarkerFaceColor = colors(3,:);
    end
    ylim([0 100]);
    xlim([0 size(num_grp,1)+1]);
    ylim([10 100]);
%     trs = find(num_grp(:,2)<=7);
%     tre = setdiff(1:size(num_grp,1),trs);
%     titl = repelem("TR",1,length(trs))';
%     titl = append(titl,string(trs));
%     titl = [titl ; append(repelem("TE",length(tre),1),string(tre'))];
    xticks(1:size(num_grp,1));
%     xticklabels(titl);
%     xlabel("Bin");
%     ylabel("percentage of cells (%)");
    box off
    % check corr
    [r1,p1] = corr((1:size(perc_mat,1))',perc_mat(:,3));
%     title({["cells grouped by activity during training"];...
%         ['Corr None: ' num2str(r1) ', p-value: ' num2str(p1)]});
    none_pnts_fit = polyfit(1:size(perc_mat,1),perc_mat(:,3),1);
    none_line = polyval(none_pnts_fit,1:size(perc_mat,1));
    plot(1:length(none_line),none_line,'Color',colors(3,:),'LineStyle','--');
    legend(ss1,"None");
end
