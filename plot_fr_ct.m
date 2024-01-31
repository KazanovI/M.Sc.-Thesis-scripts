function plot_fr_ct(FR,sess_info,ind_cells)
    global num_grp
    % binned
    fr_bin = cell(1,length(num_grp));
    [m_bin,std_bin] = deal(zeros(length(num_grp),3));
    for b = 1:length(num_grp)
        fr_bin{b} = cell(1,3);
        for s_i = 1:length(FR)
            if ismember(double(sess_info.Session(s_i)),num_grp(b,1):num_grp(b,2))
                for i = 1:size(ind_cells{s_i},2)
                    fr_bin{b}{i} = [fr_bin{b}{i}; FR{s_i}(logical(ind_cells{s_i}(:,i)))];
                end
            end
        end
        m_bin(b,:) = cellfun(@mean ,fr_bin{b});
        std_bin(b,:) = cellfun(@(x) std(x)/sqrt(length(x)),fr_bin{b});
    end
    % plot+fill per bin
    clrs = ["r","b","g"];
    figure();
    for i = 1:size(m_bin,2)
        p(i) = plot(m_bin(:,i),clrs(i));
        hold on
        x2 = [1:length(m_bin(:,i)), fliplr(1:length(m_bin(:,i)))];
        inBetween = [(m_bin(:,i)-std_bin(:,i))' fliplr((m_bin(:,i)+std_bin(:,i))')];
        fill(x2, inBetween,clrs(i),'FaceAlpha',0.3,'EdgeAlpha',0.5);
    end
    title("Mean FR per bin and cell type, all activity");
    xlabel("Bin");
    ylabel("Mean FR(Hz)");
    axis padded
    legend(p,["Enhanced","Suppressed","None"]);
    box off
    % statistics ANOVA2
    figure();
    [all_dat,cell_id,bin_id] = deal([]);
    for b = 1:length(fr_bin)
        all_dat = [all_dat ; cell2mat(fr_bin{b}')];
        cell_id = [cell_id; repelem([1 2 3]',cellfun(@length,fr_bin{b}),1)];
        bin_id = [bin_id; repelem(b,sum(cellfun(@length,fr_bin{b})),1)];
    end
    [pv,~,stats] = anovan(all_dat,{cell_id,bin_id},'model','interaction','varnames',{'Cell type','Bin'});
    [c,~,~,nms] = multcompare(stats,'Dimension',[1 2]);;
    % histograms
    figure();
    for b = 1:length(fr_bin)
        ax(b) = subplot(2,3,b);
        for clt = 1:length(fr_bin{b})
            xtk = [min(fr_bin{b}{clt}):5:max(fr_bin{b}{clt})];
            perc = (histcounts(fr_bin{b}{clt},length(xtk))/length(fr_bin{b}{clt}))*100;
            br(clt) = bar(perc,'FaceColor',clrs(clt),'FaceAlpha',0.5);
            bmean(clt) = mean(fr_bin{b}{clt});
            hold on
        end
        title(['Bin: ' num2str(b), 'Mean: ' num2str(bmean)]);
        xlabel('FR(Hz)');
        ylabel('% of cells');
        xticks(1:2:length(xtk));
        xticklabels(xtk(1:2:end));
    end
    legend(br,["Enhanced","Suppressed","None"]);
    linkaxes(ax,'xy');
end