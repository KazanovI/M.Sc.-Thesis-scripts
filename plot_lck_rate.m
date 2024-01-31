function plot_lck_rate(licks_count,trial_track)
    cnds = [1 2;3 4;5 6];
    cond_data = cell(1,length(cnds));
    for s_i = 1:length(licks_count)
        if length(unique(trial_track{s_i})) > 2
            for c_i = 1:length(cond_data)
                cond_data{c_i} = [cond_data{c_i}; mean(licks_count{s_i}(ismember(trial_track{s_i},cnds(c_i,:)),:))];
            end
        end
    end
    max_cond = cellfun(@(x) max(x(:,5:end),[],2),cond_data,'UniformOutput',false);
    m_cond = cellfun(@(x) mean(x),max_cond);
    std_cond = cellfun(@(x) std(x)/sqrt(length(x)),max_cond);
    clrs = get_color("Hit","Miss","FA");
    figure();
    b = bar(m_cond);
    b.FaceColor = 'flat';
    b.CData = clrs;
    xticklabels(["Hit","Miss","FA"]);
    hold on
    errorbar(m_cond,std_cond,'k','LineStyle','none');
    box off
    title('maximum lick probability, over frames');
    s = scatter([1 2 3],cell2mat(max_cond),20,'k','filled', 'jitter','on', 'jitterAmount',0.2);
    
    [~,~,stats] = anova1(cat(1,max_cond{:}),repelem([1 2 3]',length(max_cond{1})));
    figure();
    multcompare(stats)
end