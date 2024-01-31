function plot_start_task_svm(pa,pa_shuf,sess_info,sv)
    global num_grp
    % super-mice
    [all_res,all_res_shuf] = deal(cell(1,length(num_grp)));
    for b = 1:length(num_grp)
        for si = 1:length(pa)
            if ismember(double(sess_info.Session(si)),num_grp(b,1):num_grp(b,2))
                all_res{1,b} = [all_res{1,b} , pa{si}];
                all_res_shuf{1,b} = [all_res_shuf{1,b}, pa_shuf{si}];
            end
        end
    end
    %plot
    figure(); hold all;
    for b = 1:length(all_res)
        scatter(b,all_res{1,b},'k','filled','MarkerFaceAlpha',0.4);
        errorbar(b+0.2,mean(all_res{1,b},'omitnan'),std(all_res{1,b},[],'omitnan')/sqrt(sum(~isnan(all_res{1,b}))),'xr');
        % shuffled
        scatter(b,all_res_shuf{1,b},'b','filled','MarkerFaceAlpha',0.4);
        errorbar(b+0.2,mean(all_res_shuf{1,b},'omitnan'),std(all_res_shuf{1,b},[],'omitnan')/sqrt(sum(~isnan(all_res_shuf{1,b}))),'xr');
    end
    yl = yline(50,'-k');
    xlim([0 7]);
    ylim([20 100]);
    xticks(1:6);
    xlabel("Bin");
    ylabel("Prediction Accuracy");
    title('SVM for prediction of lick in start or task');
%     yline(mean(all_res{1,1},'omitnan'),':k')
%     yline(mean(all_res{1,2},'omitnan'),':k')
    if sv
        save2 = ['C:\onedriveitay\Master\Resnik lab\results\SVM_figs'];
        if ~exist(save2, 'dir')
            mkdir(save2);
        end
        fig_name = append("start_task",".fig");
        saveas(gcf,fullfile(save2,fig_name));
        fig_name = append("start_task",".jpeg");
        saveas(gcf,fullfile(save2,fig_name));
    end
    % anova2
    figure();
    both_res = [all_res;all_res_shuf];
    all_res_clean = cellfun(@(x) x(~isnan(x)),both_res,'UniformOutput',false);
    bin_count = cellfun(@(x) length(x), all_res_clean);
    group = [];
    for i = 1:length(bin_count)
        group = [group, repelem(i,bin_count(1,i))];
    end
    group = repmat(group,1,2);
    type = repelem([1 2],length(group)/2);
    both_dat = cell2mat(all_res_clean)';
    both_dat = both_dat(:);
    [p,tbl,stats] = anovan(both_dat,{group,type},'model','interaction','varnames',{'Bin','order'});
    figure(); 
    [c,~,~,gnms] = multcompare(stats,'Dimension',[1 2],'CType','bonferroni');
end