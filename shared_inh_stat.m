function [shrd] = shared_inh_stat(cell_inh,ind_cell_stats_task)
    inds = [cell_inh; ind_cell_stats_task];
    shrd = cell(1,length(inds));
    for s_i = 1:length(inds)
        shrd{s_i} = inds{1,s_i} == inds{2,s_i} & inds{1,s_i} == 1;
    end
    any_inh = cellfun(@(x) any(x),cell_inh);
    shrd_inh = shrd(any_inh);
    perc_sess = cellfun(@(x) (sum(x)/length(x))*100,shrd_inh,'UniformOutput',false);
%     num_sess = cellfun(@(x) sum(x),shrd_inh,'UniformOutput',false);
    dat = cell2mat(perc_sess');
    m_dat = mean(dat);
    std_dat = std(dat)/sqrt(length(dat));
    clrs = get_color(["Enhanced","Suppressed","None"]);
    figure();
    b = bar(m_dat);
    b.FaceColor = 'flat';
    b.CData = clrs;
    hold on
    errorbar(m_dat,std_dat,'k','LineStyle','none');
    ylim([0 20]);
    xticklabels(["Enhanced","Suppressed","None"]);
    ylabel("% of cells");
    title("Inhibitory cells");
    box off
end