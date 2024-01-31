function [euc_dists_it] = plot_traj(data,sess_info,inds,trial_track,phase,cond_comp,sv)
    global num_grp
    % PCA
    % by cell type, all conditions
    wind_dist = [25 29;30 34]; % calc distance between trajectories in these frames
    roi = [10 25 ; 25 40 ; 40 45 ]; % for traj plots
    tit = ["Enhanced","Suppressed","None"];
    conds_name = ["Hit","Miss","FA"];
    num_iter = 500;
    [euc_dists_it,euc_dists_it_ts] = deal(cell(1,length(cond_comp)));
    for trck = 1:length(trial_track)
        if sess_info.Phase(trck) == "Test"
            trial_track{trck}(ismember(trial_track{trck},[1 2]))  = 1;
            trial_track{trck}(ismember(trial_track{trck},[3 4]))  = 2;
            trial_track{trck}(ismember(trial_track{trck},[5 6]))  = 3;
        end
    end
    [evals_it,comps_ts_tr,comps_ts_te] = deal(cell(1,length(num_grp)));
    [evals_it{1},comps_ts_tr{1},comps_ts_te{1}] = deal(cell(1,size(inds{1},2)));
    [var_exp,max_var_exp,var_exp_perc] = deal(zeros(length(num_grp),size(inds{1},2),num_iter));
    num_cells = zeros(length(num_grp),size(inds{1},2));
    for it = 1:num_iter
        for ct = 1:size(inds{1},2) % cell type
            euc_dists = cell(1,length(num_grp)); % euclidian distance in test
            for b = 1:size(num_grp,1) % per bin
                all_cells = [];
                all_cells_cond = cell(1,2);
                for s_i = 1:length(data)% session
                    if ismember(double(sess_info.Session(s_i)),num_grp(b,1):num_grp(b,2))
                        dat_cell = data{s_i}(logical(inds{s_i}(:,ct)),:,:);
                        if ~isempty(dat_cell) & ~all(isnan(dat_cell),'all')
                            dat_cell(std(dat_cell,[],[2 3]) > 3,:,:) = [];
                            if sess_info.Phase(s_i) == "Test"
                                trialsN = histcounts(trial_track{s_i},1:1:4);
                                if ~any(trialsN(cond_comp) < 5) | size(dat_cell,1) == 0
                                    ntrials_cond_comp = trialsN(cond_comp);
                                    [~,m_ind] = min(ntrials_cond_comp); % condition with less trials
                                    tr_ind1 = randsample(ntrials_cond_comp(m_ind),ntrials_cond_comp(m_ind),true); % randomly take trials
                                    tr_ind2 = randsample(ntrials_cond_comp(3-m_ind),ntrials_cond_comp(m_ind),true); % randomly take trials
                                    tr_cond_ind = find(trial_track{s_i} == cond_comp(m_ind)); % check indices and extract from each cond
                                    tr_cond_ind = tr_cond_ind(tr_ind1); % smaller group
                                    tr1 = dat_cell(:,:,tr_cond_ind);
                                    tr_cond_ind2 = find(trial_track{s_i} == cond_comp(3-m_ind));
                                    tr_cond_ind2 = tr_cond_ind2(tr_ind2);
                                    tr2 = dat_cell(:,:,tr_cond_ind2);
                                    t_cat = [cat(3,tr1, tr2)];
                                    for cond = 1:length(cond_comp)
                                        if cond == m_ind
                                            mean_cond = squeeze(mean(tr1,3)); % mean over trials
                                        else
                                            mean_cond = squeeze(mean(tr2,3)); % mean over trials
                                        end
                                        if min(size(mean_cond)) == 1 & ~isrow(mean_cond)
                                            mean_cond = mean_cond';
                                        end
                                        all_cells_cond{cond} = [all_cells_cond{cond}; mean_cond]; % concat
                                    end
                                else
                                    continue;
                                end
                            else
                                tr_ind = randsample(size(dat_cell,3),size(dat_cell,3),true); 
                                t_cat = dat_cell(:,:,tr_ind);
                            end
                            t1 = squeeze(mean(t_cat,3)); % mean over trials
                            if ~all(isnan(t1),'all')
                                all_cells = [all_cells ; t1];
                            end
                        end
                    end
                end
                if it == 1
                    num_cells(b,ct) = size(all_cells,1);
                end
                m_all = mean(all_cells,2);
                all_cells_m = all_cells - m_all; % mean center
                roi_cov = cov(all_cells_m'); % cov mat by all data
                % eig and sort
                [evecs,evals] = eig( roi_cov );
                [evals,sidx] = sort(diag(evals),'descend');
                evecs = evecs(:,sidx);
                evals = 100*(evals/sum(evals));
                % var data
                num_pcs = find(cumsum(evals) >= 80);
                var_exp(b,ct,it) = num_pcs(1); % how many components to explain 80%
                var_exp_perc(b,ct,it) = (num_pcs(1)/length(evals))*100;
                max_var_exp(b,ct,it) = evals(1);
                % principal components time series first 3 otrhogonal comps
                if b > 3
                    % map to data
                    dat_hit = all_cells_cond{1};
                    dat_miss = all_cells_cond{2};
                    dat_hit = dat_hit - m_all;
                    dat_miss = dat_miss - m_all;
                    pc_timeseries_hit = evecs(:,1:3)'*dat_hit;
                    pc_timeseries_miss = evecs(:,1:3)'*dat_miss;
                    dat_c = [pc_timeseries_hit' pc_timeseries_miss'];

                    % for dist measures
                    if var_exp(b,ct,it) > 3
                        h = evecs(:,1:var_exp(b,ct,it))'*dat_hit;
                        m = evecs(:,1:var_exp(b,ct,it))'*dat_miss;
                    else
                        h = evecs(:,1:3)'*dat_hit;
                        m = evecs(:,1:3)'*dat_miss; 
                    end
                    pcs_for_dist = cat(3,h,m);
                else
                    dat_c = evecs(:,1:3)'*all_cells_m;
                end
                if b > 3
                    pc_conds = cat(3,pc_timeseries_hit,pc_timeseries_miss);
                else
                    pc_conds = dat_c;
                end
                pc_timeseries = smoothdata(pc_conds,2,'movmean',10); % smooth
                % plots
                if it == 1
                    if b > 3
                        comps_ts_te{b}{ct} = dat_c;
                    else
                        comps_ts_tr{b}{ct} = dat_c;
                    end
                    try
                        evals_it{b}{ct} = evals(1:20);
                    catch
                        evals_it{b}{ct} = evals(1:end);
                    end
                elseif it ~= num_iter
                    if b > 3
                        comps_ts_te{b}{ct} = comps_ts_te{b}{ct} + dat_c;
                    else
                        comps_ts_tr{b}{ct} = comps_ts_tr{b}{ct} + dat_c;
                    end
                    evals_it{b}{ct} = evals_it{b}{ct} + evals(1:length(evals_it{b}{ct}));
                else
                    if b > 3
                        comps_ts_te{b}{ct} = comps_ts_te{b}{ct}./num_iter;
                    else
                        comps_ts_tr{b}{ct} = comps_ts_tr{b}{ct}./num_iter;
                    end
                    evals_it{b}{ct} = evals_it{b}{ct}./num_iter;
                    ln = ["-","--","-."];
                    
                    figure();
                    subplot(2,1,1)
                    if b > 3
                        fp = plot(1:size(comps_ts_te{b}{ct},1),comps_ts_te{b}{ct}');
                        xticks(0:5:size(comps_ts_te{b}{ct},1));
                    else
                        fp = plot(1:size(comps_ts_tr{b}{ct},2),comps_ts_tr{b}{ct}');
                        xticks(0:5:size(comps_ts_tr{b}{ct},2));
                    end
%                     title("PC activity");
                    ylabel("magnitude");
                    xlabel("Time relative to lick onset(Sec)");
                    digits(2);
                    xticklabels(string(vpa(-30:5:15)/30));
                    if b > 3
                        leg = append(conds_name(cond_comp),["1","2","3"]');
                        leg = leg(:);
                        legend(leg,'Location','northwest');
                        fp(1).LineStyle = ln(1); fp(4).LineStyle = ln(1);
                        fp(2).LineStyle = ln(2); fp(5).LineStyle = ln(2);
                        fp(3).LineStyle = ln(3); fp(6).LineStyle = ln(3);
                    else
                        legend("PC1","PC2","PC3",'Location','northwest');
                        fp(1).LineStyle = ln(1);
                        fp(2).LineStyle = ln(2);
                        fp(3).LineStyle = ln(3);
                    end
                    xline(30,'k','LineWidth',2,'DisplayName','Lick time');
                    box off
                    subplot(2,1,2)
                    plot(cumsum(evals_it{b}{ct}),'ko');
                    ylim([0 100]);
                    xlim([0 length(evals_it{b}{ct})]);
                    yline(80,':k');
%                     title("Scree plot");
                    xlabel("Components");
                    ylabel("Cumulative variance explained(%)");
%                     sgtitle([char(phase) ', Bin: ' num2str(b)]);
                    box off
                    % save figs
                    if sv
                        save2 = ['C:\onedriveitay\Master\Resnik lab\results\PCA figs\celltype\' char(phase)];
                        if ~exist(save2, 'dir')
                            mkdir(save2);
                        end
                        fig_name = append("Comps_",tit(ct),num2str(b),"_",conds_name(cond_comp(1)),"_",conds_name(cond_comp(2)),".fig");
                        saveas(gcf,fullfile(save2,fig_name));
                        fig_name = append("Comps_",tit(ct),num2str(b),"_",conds_name(cond_comp(1)),"_",conds_name(cond_comp(2)),".jpeg");
                        saveas(gcf,fullfile(save2,fig_name));
                    end
                    close;
                    
                    % state space
                    clrs = get_color(conds_name(cond_comp));
                    figure(); hold on;
                    for reg = 1:size(roi,1) % temporal regions
                        for cond = 1:size(pc_timeseries,3)
                            cond_state = pc_timeseries(:,roi(reg,1):roi(reg,2),cond);
                            if reg ~= 2
                                plot3(cond_state(1,:),cond_state(2,:),cond_state(3,:),'Color',clrs(cond,:),'LineStyle',':','LineWidth',1);
                                if reg == 1
                                    plot3(cond_state(1,1),cond_state(2,1),cond_state(3,1),'ko');
                                elseif reg == 3
                                    plot3(cond_state(1,end),cond_state(2,end),cond_state(3,end),'kp');
                                end
                            else
                                plot3(cond_state(1,:),cond_state(2,:),cond_state(3,:),'Color',clrs(cond,:),'LineStyle','-','LineWidth',2);
                                plot3(cond_state(1,1),cond_state(2,1),cond_state(3,1),'ks');
                            end
                        end
                        xlabel('PC1');
                        ylabel('PC2');
                        zlabel('PC3');
                        if b > 3
                            legend(conds_name(cond_comp(1)),"",conds_name(cond_comp(2)),"",'Location','best');
                        end
                        view(3); 
                        rotate3d on
                        grid on
                    end
                    % save figs
                    if sv
                        save2 = ['C:\onedriveitay\Master\Resnik lab\results\PCA figs\celltype\' char(phase)];
                        if ~exist(save2, 'dir')
                            mkdir(save2);
                        end
                        fig_name = append(tit(ct),num2str(b),"_",conds_name(cond_comp(1)),"_",conds_name(cond_comp(2)),".fig");
                        saveas(gcf,fullfile(save2,fig_name));
                        fig_name = append(tit(ct),num2str(b),"_",conds_name(cond_comp(1)),"_",conds_name(cond_comp(2)),".jpeg");
                        saveas(gcf,fullfile(save2,fig_name));
                    end
                    close all;
                end
                if b > 3
                    % Euclidian distance
                    d = zeros(1,size(pcs_for_dist,2));
                    pnts = zeros(size(pcs_for_dist,[1 2]));
                    for p = 1:size(pcs_for_dist,1) % num pc's
                        pnts(p,:) = (pcs_for_dist(p,:,1) - pcs_for_dist(p,:,2)).^2;
                    end
                    d = sqrt(sum(pnts));
                    euc_dists{b}{1} = d;
                    euc_dists{b}{2} = [mean(d(wind_dist(1,1):wind_dist(1,2))) mean(d(wind_dist(2,1):wind_dist(2,2)))];
                    euc_dists{b}{3} = mean(d(wind_dist(2,1):wind_dist(2,2))) - mean(d(wind_dist(1,1):wind_dist(1,2)));
                end
            end
            close all
            euc_dists_it{ct}(it,:) = cellfun(@(x) x{3},euc_dists(4:end));
            euc_dists_it_ts{ct}(it,:) = cell2mat(cellfun(@(x) x{1},euc_dists(4:end),'UniformOutput',false));
        end
    end
    % plot results
    mean_dists = cell2mat(cellfun(@(x) mean(x),euc_dists_it,'UniformOutput',false)');
    ste_dists = cell2mat(cellfun(@(x) std(x)/sqrt(length(x)),euc_dists_it,'UniformOutput',false)');
    mean_ts = cellfun(@(x) reshape(x,[],3)',cellfun(@(x) mean(x),euc_dists_it_ts,'UniformOutput',false),'UniformOutput',false);
    ste_ts = cellfun(@(x) reshape(x,[],3)',cellfun(@(x) std(x),euc_dists_it_ts,'UniformOutput',false),'UniformOutput',false);
    if phase == "mix"
        clrs = get_color("mix Enhanced","mix Suppressed","mix None");
    else
        clrs = get_color("Enhanced","Suppressed","None");
    end
    ofst = [-0.2 0 0.2];
    figure();
    for c = 1:length(mean_dists)
        bar([1 2 3] + ofst(c),mean_dists(c,:),'FaceColor',clrs(c,:),'BarWidth',0.2);
        hold on
        errorbar([1 2 3] + ofst(c),mean_dists(c,:),ste_dists(c,:),'LineStyle','none','Color',[.5 .5 .5]);
    end
    ylabel('Diff distance');
    xlabel('Bin');
    xticks(1:3)
    xticklabels(["4","5","6"]);
%     ylim([-1 2]);
    box off
    if sv
        % save figs
        fig_name = append("barDists_",conds_name(cond_comp(1)),"_",conds_name(cond_comp(2)),".fig");
        saveas(gcf,fullfile(save2,fig_name));
        fig_name = append("barDists_",conds_name(cond_comp(1)),"_",conds_name(cond_comp(2)),".jpeg");
        saveas(gcf,fullfile(save2,fig_name));
    end
    close;
    
    figure();
    ln = ["-","--","-."];
    clear fp
    for c = 1:length(mean_ts)
        vecs = mean_ts{c};
        vecs_sd = ste_ts{c};
%         dists = mean_dists(c,:);
        ax(c) = subplot(3,1,c);
        for bin = 1:3
            vec_bin = vecs(bin,:);
            vec_sd_bin = vecs_sd(bin,:);
            fp(bin) = plot(1:length(vec_bin),vec_bin,'Color',clrs(c,:),'LineStyle',ln(bin));
            hold on 
            x2 = [1:length(vec_bin), fliplr(1:length(vec_bin))];
            inBetween = [(vec_bin-vec_sd_bin) fliplr((vec_bin+vec_sd_bin))];
            fill(x2, inBetween,clrs(c,:),'FaceAlpha',0.3,'EdgeAlpha',0.5);
        end
        box off
        ylabel("Distance");
        xticks(0:5:length(vec_bin));
        digits(2);
        xticklabels(string(vpa(-30:5:15)/30));
%         title({tit(c),num2str(dists)});
        if c == 1
            legend(fp,"TE4","TE5","TE6",'Location','southwest');
        end
    end
    xlabel('Time relative to lick onset(s)');
    linkaxes(ax,'y');
%     cnds = conds_name(cond_comp);
%     cnds_name = append(cnds(1)," and ",cnds(2));
%     sgtitle([char(phase) ', dist of ' char(cnds_name) ',Num iters: ' num2str(num_iter)]);
    if sv
        % save figs
        fig_name = append("ts_Bins_",conds_name(cond_comp(1)),"_",conds_name(cond_comp(2)),".fig");
        saveas(gcf,fullfile(save2,fig_name));
        fig_name = append("ts_Bins_",conds_name(cond_comp(1)),"_",conds_name(cond_comp(2)),".jpeg");
        saveas(gcf,fullfile(save2,fig_name));
    end
    close;
    % same but mean over bins
    figure();
    clear fp
    for c = [3 2 1]
        vecs = mean_ts{c};
%         dists = mean_dists(c,:);
%         ax(c) = subplot(3,1,c);
        vec_bin = mean(vecs);
        vec_sd_bin = std(vecs)/sqrt(size(vecs,1));
        plot(1:length(vec_bin),vec_bin,'Color',clrs(c,:),'LineWidth',3);
        hold on 
        x2 = [1:length(vec_bin), fliplr(1:length(vec_bin))];
        inBetween = [(vec_bin-vec_sd_bin) fliplr((vec_bin+vec_sd_bin))];
        fill(x2, inBetween,clrs(c,:),'FaceAlpha',0.3,'EdgeAlpha',0.3);

    end
        box off
        ylabel("Distance");
        xticks(0:5:length(vec_bin));
        digits(2);
        xticklabels(string(vpa(-30:5:15)/30));
        xlim([25 35]); 
        ylim([2 10]);
    xlabel('Time relative to lick onset(s)');
%     linkaxes(ax,'y');
%     cnds = conds_name(cond_comp);
%     cnds_name = append(cnds(1)," and ",cnds(2));
%     sgtitle([char(phase) ', dist of ' char(cnds_name) ',Num iters: ' num2str(num_iter)]);
    if sv
        % save figs
        fig_name = append("ts_Bins_m",conds_name(cond_comp(1)),"_",conds_name(cond_comp(2)),".fig");
        saveas(gcf,fullfile(save2,fig_name));
        fig_name = append("ts_Bins_m",conds_name(cond_comp(1)),"_",conds_name(cond_comp(2)),".jpeg");
        saveas(gcf,fullfile(save2,fig_name));
    end
    close;
    % 
%     for bn = 1:length(euc_dists_it)
%         dat_bin = cell2mat(cellfun(@(x) x(:,bn),euc_dists_it,'UniformOutput',false));
%         figure();
%         for dt = 1:size(dat_bin,2)-1
%             histogram(dat_bin(:,dt),'FaceColor',clrs(dt,:),'Normalization','probability');
%             hold on;
%         end
%         xlabel("Distance(After-Before)");
%         ylabel("Probabilty");
%         legend(tit);
% %         title(['Bin: ' num2str(bn)]);
%         box off
%         if sv
%             % save figs
%             fig_name = append("dists_",num2str(bn),"_",conds_name(cond_comp(1)),"_",conds_name(cond_comp(2)),".fig");
%             saveas(gcf,fullfile(save2,fig_name));
%             fig_name = append("dists_",num2str(bn),"_",conds_name(cond_comp(1)),"_",conds_name(cond_comp(2)),".jpeg");
%             saveas(gcf,fullfile(save2,fig_name));
%         end
%         close;
%     end

%     % look at the distribution of dists
%     for i = 1:length(euc_dists_it)
%         figure(); hold on;
%         for bn = 1:3
%             histogram(euc_dists_it{i}(:,bn));
%         end
%         legend("TE4","TE5","TE6",'Location','northeast');
%         title(tit(i));
%     end
%     [~,p] = ttest2(euc_dists_it{2}(:,2),euc_dists_it{2}(:,3));
    % confidence intervals, effect of time
%     figure();
%     histogram(euc_dists_it{1}(:,3) - euc_dists_it{1}(:,1)); hold on;
%     perc = prctile(euc_dists_it{1}(:,3) - euc_dists_it{1}(:,1),[2.5 97.5]);
%     histogram(euc_dists_it{2}(:,3) - euc_dists_it{2}(:,1));
%     histogram(euc_dists_it{3}(:,3) - euc_dists_it{3}(:,1));
%     title("TE6-TE4 diff of distances of all cell types")
%     xline(perc,':r','LineWidth',2);
%     legend([tit,"2.5%","97.5%"]);
    % [~,p] = kstest2(euc_dists_it{2}(:,3) - euc_dists_it{2}(:,1),euc_dists_it{3}(:,3) - euc_dists_it{3}(:,1));

%     % confidence intervals, effect of celltype,bin6
%     figure();
%     h1 = histogram(euc_dists_it{1}(:,3) - euc_dists_it{2}(:,3),'Normalization','probability'); hold on;
%     perc = prctile(euc_dists_it{1}(:,3) - euc_dists_it{2}(:,3),[2.5 97.5]);
%     h2 = histogram(euc_dists_it{1}(:,3) - euc_dists_it{3}(:,3),'Normalization','probability');
%     title("TE6-TE6 diff of distances of all cell types")
%     x1 = xline(perc,':r','LineWidth',2);
%     legend([h1; h2; x1],[tit(1:2),"CI"]);
end