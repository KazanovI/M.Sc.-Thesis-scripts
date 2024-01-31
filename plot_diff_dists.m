function plot_diff_dists(dists_task,dists_mix,cond_comp)
    conds_name = ["Hit","Miss","FA"];
    both = [dists_mix; dists_task];
    pnts = cell(1,size(both,1));
    for i = 1:size(both,1)
        mean_dists(i,:) = mean(cell2mat(cellfun(@(x) mean(x),both(i,:),'UniformOutput',false)'),2);
        ste_dists(i,:) = mean(cell2mat(cellfun(@(x) std(x)/sqrt(length(x)),both(i,:),'UniformOutput',false)'),2);
        pnts{i} = cell2mat(cellfun(@(x) mean(x,2),both(i,:),'UniformOutput',false));
    end
    
    clrs_task = get_color("Enhanced","Suppressed","None");
    clrs_mix = get_color("mix Enhanced","mix Suppressed","mix None");
    figure();
    b = bar([1 2 3],mean_dists);
    b(1).FaceColor = 'flat';
    b(1).CData = clrs_mix;
    b(2).FaceColor = 'flat';
    b(2).CData = clrs_task;
    hold on
    for b_i = 1:length(b)
        xend = b(b_i).XEndPoints;
        errorbar(xend,mean_dists(b_i,:),ste_dists(b_i,:),'k','LineStyle','none');
%         scatter(xend,pnts{b_i},15,'k','filled');
    end
    xticklabels(["Enhanced","Suppressed","None"]);
    box off
    % save figs
    save2 = ['C:\onedriveitay\Master\Resnik lab\results\PCA figs\celltype'];
    fig_name = append("barDists_both",conds_name(cond_comp(1)),"_",conds_name(cond_comp(2)),".fig");
    saveas(gcf,fullfile(save2,fig_name));
    fig_name = append("barDists_both",conds_name(cond_comp(1)),"_",conds_name(cond_comp(2)),".jpeg");
    saveas(gcf,fullfile(save2,fig_name));
    % anova2 cell type and time
    data_all = cellfun(@(x) mean(x,2),both,'UniformOutput',false);
    ct = repmat(repelem([1 2 3]',500),2,1);
    time = repelem([1 2]',1500);
    data_all = data_all';
    [~,~,stats] = anovan(cat(1,data_all{:}),{ct time},'model',2,'varnames',{'cell type','time'});
    figure();
    [c,~,~,gnms] = multcompare(stats,'Dimension',[1 2])
    % permutations
    dat_mob = cellfun(@(x) mean(x,2),both,'UniformOutput',false); % mean over bins
    num_iters = 5000;
    ct_names = ["Enhanced","Suppressed","None"];
    for i = 1:length(dat_mob) % per cell type
        pv = permutationTest(dat_mob{1,i}, dat_mob{2,i},num_iters,'sidedness','smaller');
        fprintf('\n %s ,dist: mix vs. task. p-value: %.4f, iters: %d \n',ct_names(i),pv,num_iters);
    end
    % by bin 
    for ti = 1:2
        for bn = 1:length(both)
            dat_bins = cell2mat(cellfun(@(x) x(:,bn),both(ti,1:2),'UniformOutput',false));
            figure();
            if ti == 1
                clrs = clrs_mix;
            else
                clrs = clrs_task;
            end
            for dt = 1:size(dat_bins,2)
                histogram(dat_bins(:,dt),'FaceColor',clrs(dt,:),'Normalization','probability');
                hold on;
            end
            xlabel("Distance(After-Before)");
            ylabel("Probabilty");
%             legend(ct_names);
    %         title(['Bin: ' num2str(bn)]);
            box off
            pv = permutationTest(dat_bins(:,1), dat_bins(:,2),num_iters,'sidedness','larger');
            fprintf('\n dist: enh vs. sup. p-value: %.4f, iters: %d \n',pv,num_iters);
%             close;
        end
    end
end