function plot_perc_var_exp(start,task)
    global num_grp
    dat = [mean(start,2) mean(task,2)];
    m_dat = mean(dat);
    std_dat = std(dat)/sqrt(size(dat,1));
    figure();
    bar(m_dat,'w');
    hold on
    errorbar(m_dat,std_dat,'k','LineStyle','none');
    ylabel('% components');
    xticklabels(["Pre-task" "Task"]);
    title("% of components to explain 80% variance");
    box off
    [h,p] = ttest2(dat(:,1),dat(:,2));
end