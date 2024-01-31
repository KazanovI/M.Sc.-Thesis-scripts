function [var_exp_perc] = plot_traj_all(data,sess_info,trial_track,phase,cond_comp,sv)
    global num_grp
    %% all cells ,hit&miss, all imaged cells concatenated by bin
    % to see learning manifested somehow through bins
    num_iter = 500;
    wind_dist = [25 29;30 34]; % calc distance between trajectories in these frames
    roi = [10 25 ; 25 40 ; 40 45 ]; % for traj plots
    conds_name = ["Hit","Miss","FA"];
    euc_dists_it = zeros(num_iter,3); % distance
    [euc_dists2_it] = deal(zeros(num_iter,length(num_grp))); 
    euc_dists_it_ts = zeros(num_iter,135); % time-series
    for trck = 1:length(trial_track)
        if sess_info.Phase(trck) == "Test"
            trial_track{trck}(ismember(trial_track{trck},[1 2]))  = 1;
            trial_track{trck}(ismember(trial_track{trck},[3 4]))  = 2;
            trial_track{trck}(ismember(trial_track{trck},[5 6]))  = 3;
        end
    end
    pc_it = cell(1,length(num_grp));
    [evals_it,comps_ts_tr,comps_ts_te] = deal(cell(1,length(num_grp)));
    [var_exp,var_exp_perc,max_var_exp,var_exp_perc2] = deal(zeros(length(num_grp),num_iter));
    num_cells = zeros(1,length(num_grp)); % how many cells in bin
    for it = 1:num_iter
        euc_dists = cell(1,length(num_grp)); % euclidian distance in test
        for b = 1:size(num_grp,1) % per bin
            all_cells = [];
            all_cells_cond = cell(1,2);
            for s_i = 1:length(data)% session
                if ismember(double(sess_info.Session(s_i)),num_grp(b,1):num_grp(b,2))
                    dat_sess = data{s_i};
                    dat_sess(std(dat_sess,[],[2 3]) > 3,:,:) = [];
                    if sess_info.Phase(s_i) == "Test"
                        trialsN = histcounts(trial_track{s_i},1:1:4);
                        if ~any(trialsN(cond_comp) < 5)
                            ntrials_cond_comp = trialsN(cond_comp);
                            [~,m_ind] = min(ntrials_cond_comp); % condition with less trials
                            tr_ind1 = randsample(ntrials_cond_comp(m_ind),ntrials_cond_comp(m_ind),true); % randomly take trials
                            tr_ind2 = randsample(ntrials_cond_comp(3-m_ind),ntrials_cond_comp(m_ind),true); % randomly take trials
                            tr_cond_ind = find(trial_track{s_i} == cond_comp(m_ind)); % check indices and extract from each cond
                            tr_cond_ind = tr_cond_ind(tr_ind1); % smaller group
                            tr1 = dat_sess(:,:,tr_cond_ind);
                            tr_cond_ind2 = find(trial_track{s_i} == cond_comp(3-m_ind));
                            tr_cond_ind2 = tr_cond_ind2(tr_ind2);
                            tr2 = dat_sess(:,:,tr_cond_ind2);
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
                        tr_ind = randsample(size(dat_sess,3),size(dat_sess,3),true); 
                        t_cat = dat_sess(:,:,tr_ind); % train
                    end
                    t1 = squeeze(mean(t_cat,3)); % mean over trials
                    all_cells = [all_cells ; t1];
                end
            end
            if it == 1
                num_cells(b) = size(all_cells,1);
            end
            m_val_cells = mean(all_cells,2);
            all_cells_m = all_cells - m_val_cells; % mean center
            roi_cov = cov(all_cells_m'); % cov mat by all data

            % eig and sort
            [evecs,evals] = eig( roi_cov );
            [evals,sidx] = sort(diag(evals),'descend');
            evecs = evecs(:,sidx);
            evals = 100*evals/sum(evals);
            % var data
            num_pcs = find(cumsum(evals) >= 80);
            var_exp(b,it) = num_pcs(1);
            var_exp_perc(b,it) = (num_pcs(1)/length(evals))*100;
            var_exp_perc2(b,it) = sum(evals(1:3));
            max_var_exp(b,it) = evals(1);
            % principal components time series first 3 otrhogonal comps
            if b > 3
                % map to data
                dat_hit = all_cells_cond{1};
                dat_miss = all_cells_cond{2};
                dat_hit = dat_hit - m_val_cells;
                dat_miss = dat_miss - m_val_cells;
                pc_timeseries_hit = evecs(:,1:3)'*dat_hit;
                pc_timeseries_miss = evecs(:,1:3)'*dat_miss;
                dat_c = [pc_timeseries_hit' pc_timeseries_miss'];

                % for dist measures
                h = evecs(:,1:var_exp(b))'*dat_hit;
                m = evecs(:,1:var_exp(b))'*dat_miss;
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

%             plots
            if it == 1
                if b > 3
                    comps_ts_te{b} = dat_c;
                else
                    comps_ts_tr{b} = dat_c;
                end
                evals_it{b} = evals(1:20);
            elseif it ~= num_iter
                if b > 3
                    comps_ts_te{b} = comps_ts_te{b} + dat_c;
                else
                    comps_ts_tr{b} = comps_ts_tr{b} + dat_c;
                end
                evals_it{b} = evals_it{b} + evals(1:20);
            else
                comps_ts_te{b} = comps_ts_te{b}./num_iter;
                comps_ts_tr{b} = comps_ts_tr{b}./num_iter;
                evals_it{b} = evals_it{b}./num_iter;
                ln = ["-","--","-."];
                figure();
                subplot(2,1,1)
                if b > 3
                    fp = plot(1:size(comps_ts_te{b},1),comps_ts_te{b}');
                    xticks(0:5:size(comps_ts_te{b},1));
                else
                    fp = plot(1:size(comps_ts_tr{b},2),comps_ts_tr{b}');
                    xticks(0:5:size(comps_ts_tr{b},2));
                end
%                 title("PC activity");
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
                plot(cumsum(evals_it{b}),'ko'); ylim([0 100]);
%                 title("Scree plot");
                xlabel("Components");
                ylabel("Cumulative variance explained(%)");
                yline(80,':k');
%                 sgtitle(['all cells, Bin: ' num2str(b)]);
                box off
                if sv
                    save2 = ['C:\onedriveitay\Master\Resnik lab\results\PCA figs\All cells\' char(phase)];
                    if ~exist(save2, 'dir')
                        mkdir(save2);
                    end
                    fig_name = append("Comps_",num2str(b),"_",conds_name(cond_comp(1)),"_",conds_name(cond_comp(2)),".fig");
                    saveas(gcf,fullfile(save2,fig_name));
                    fig_name = append("Comps_",num2str(b),"_",conds_name(cond_comp(1)),"_",conds_name(cond_comp(2)),".jpeg");
                    saveas(gcf,fullfile(save2,fig_name));
                end
                close;
            end
            
            if it == 1  
                pc_it{b} = pc_timeseries;
            elseif it < num_iter
                pc_it{b} = pc_it{b} + pc_timeseries;
            elseif it == num_iter
                pc_it{b} = pc_it{b} + pc_timeseries;
                pc_it{b} = pc_it{b}./num_iter;
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
        euc_dists_it(it,:) = cellfun(@(x) x{3},euc_dists(4:end)); % diff distance value
        euc_dists2_it(it,:) = cell2mat(cellfun(@(x) x{2},euc_dists(4:end),'UniformOutput',false));% distance values base and lick
        euc_dists_it_ts(it,:) = cell2mat(cellfun(@(x) x{1},euc_dists(4:end),'UniformOutput',false));
    end
    %state space
    for b = 1:length(num_grp)
        if b > 3
            clrs = get_color(conds_name(cond_comp));
        else
            clrs = [0.1 0.1 0.1];
        end
        figure(); hold on;
        for reg = 1:3 % temporal regions
            for cond = 1:size(pc_it{b},3)
                cond_state = pc_it{b}(:,roi(reg,1):roi(reg,2),cond);
                if reg ~= 2
                    plot3(cond_state(1,:),cond_state(2,:),cond_state(3,:),'Color',clrs(cond,:),'LineStyle',':','LineWidth',1);
                    if reg == 1
                        plot3(cond_state(1,1),cond_state(2,1),cond_state(3,1),'ko');
                    elseif reg == 3
                        plot3(cond_state(1,end),cond_state(2,end),cond_state(3,end),'ks');
                    end
                else
                    plot3(cond_state(1,:),cond_state(2,:),cond_state(3,:),'Color',clrs(cond,:),'LineStyle','-','LineWidth',2);
%                     plot3(cond_state(1,1),cond_state(2,1),cond_state(3,1),'ks');
                end
            end
            xlabel('PC1');
            ylabel('PC2');
            zlabel('PC3');
%             title(['Bin: ' num2str(b)]);
            if b > 3
                legend(conds_name(cond_comp(1)),"",conds_name(cond_comp(2)),"");
            end  
        end
        view(3);
        rotate3d on
        grid on
        if sv
            save2 = ['C:\onedriveitay\Master\Resnik lab\results\PCA figs\All cells\' char(phase)];
            if ~exist(save2, 'dir')
                mkdir(save2);
            end
            fig_name = append(num2str(b),"_",conds_name(cond_comp(1)),"_",conds_name(cond_comp(2)),".fig");
            saveas(gcf,fullfile(save2,fig_name));
            fig_name = append(num2str(b),"_",conds_name(cond_comp(1)),"_",conds_name(cond_comp(2)),".jpeg");
            saveas(gcf,fullfile(save2,fig_name));
        end
    end
    close all
    
    % plot results
    mean_dists = mean(euc_dists_it);
    ste_dists = std(euc_dists_it)/sqrt(size(euc_dists_it,1));
    mean_ts =  reshape(mean(euc_dists_it_ts),[],3)';
    ste_ts =  reshape(std(euc_dists_it_ts)/sqrt(length(euc_dists_it_ts)),[],3)';
    mean2 =  reshape(mean(euc_dists2_it),[],3)';
    ste2 = reshape(std(euc_dists2_it)/sqrt(size(euc_dists2_it,1)),[],3)';
    
    figure();
    bar([1 2 3],mean_dists,'FaceColor','k');
    hold on
    errorbar([1 2 3],mean_dists,ste_dists,'LineStyle','none','Color',[.5 .5 .5]);
    ylabel('Diff distance');
    xlabel('Bin');
    xticklabels(["4","5","6"]);
    ylim([-1 2]);
    box off
    if sv
        % save figs
        fig_name = append("barDists_",conds_name(cond_comp(1)),"_",conds_name(cond_comp(2)),".fig");
        saveas(gcf,fullfile(save2,fig_name));
        fig_name = append("barDists_",conds_name(cond_comp(1)),"_",conds_name(cond_comp(2)),".jpeg");
        saveas(gcf,fullfile(save2,fig_name));
    end
    % before and after
    figure();
    br = bar([1 2 3],mean2,'FaceColor','k');
    hold on
    errorbar([br(1).XEndPoints' br(2).XEndPoints'],mean2,ste2,'LineStyle','none','Color',[.5 .5 .5]);
    ylabel('distance');
    xlabel('Bin');
    xticklabels(["4","5","6"]);
    ylim([4 7]);
    xpnt = [br(1).XEndPoints' br(2).XEndPoints']';
    scatter(reshape(repelem(xpnt(:),500),500,6),euc_dists2_it,12,'MarkerFaceColor','w',...
        'jitter','on','jitterAmount',0.05,'MarkerEdgeColor','k');
    ylim([3.5 11]);
    box off
    if sv
        % save figs
        fig_name = append("barDists2_",conds_name(cond_comp(1)),"_",conds_name(cond_comp(2)),".fig");
        saveas(gcf,fullfile(save2,fig_name));
        fig_name = append("barDists2_",conds_name(cond_comp(1)),"_",conds_name(cond_comp(2)),".jpeg");
        saveas(gcf,fullfile(save2,fig_name));
    end
%     grps = repelem([1 2 3]',length(euc_dists_it));
%     anovan(euc_dists_it(:),grps);
    close;
    figure();
    for i = 1:size(mean_ts,1)
        ax(i) = subplot(3,1,i);
        plot(1:length(mean_ts(i,:)),mean_ts(i,:),'k');
        hold on
        x2 = [1:length(mean_ts(i,:)), fliplr(1:length(mean_ts(i,:)))];
        inBetween = [(mean_ts(i,:)-ste_ts(i,:)) fliplr((mean_ts(i,:)+ste_ts(i,:)))];
        fill(x2, inBetween,[0.2 0.2 0.2],'FaceAlpha',0.3,'EdgeAlpha',0.5);
        xticks(0:5:size(mean_ts,2));
        digits(2);
        xticklabels(string(vpa(-30:5:15)/30));
%         title(['Dist lick(Lick-Base): ' num2str(mean_dists(i))]);
        ylabel("Distance");
        box off
    end
    linkaxes(ax,'y');
%     sgtitle({[char(phase) ', All cells, distance of ' char(conds_name(cond_comp(1))) ' and '...
%         char(conds_name(cond_comp(2))) ' trajectories'],['Num Iters: ' num2str(num_iter)]});
    xlabel("Time relative to lick onset(Sec)");
    if sv
        % save figs
        fig_name = append("ts_Bins_",conds_name(cond_comp(1)),"_",conds_name(cond_comp(2)),".fig");
        saveas(gcf,fullfile(save2,fig_name));
        fig_name = append("ts_Bins_",conds_name(cond_comp(1)),"_",conds_name(cond_comp(2)),".jpeg");
        saveas(gcf,fullfile(save2,fig_name));
    end
    close;
    figure();
    for i = 1:size(mean2,1)
        axx(i) = subplot(size(mean2,1),1,i);
        barh([0.8 1.2],mean2(i,:),'k','BarWidth',0.4);
        hold on
        errorbar(mean2(i,:),[0.8 1.2],ste2(i,:),'horizontal','LineStyle','none','Color',[.5 .5 .5]');
        yticks([0.8 1.2]);
        yticklabels({"Base","Lick"});
%         title(['Bin: ' num2str(i+3)]);
        box off
    end
    xlabel("Dist");
%     sgtitle(append(phase,", Distance between " ,conds_name(cond_comp(1)), " and " ,conds_name(cond_comp(2)), " around lick"))
    linkaxes(axx,'x');
    axis padded
    if sv
        % save figs
        fig_name = append("dists_",conds_name(cond_comp(1)),"_",conds_name(cond_comp(2)),".fig");
        saveas(gcf,fullfile(save2,fig_name));
        fig_name = append("dists_",conds_name(cond_comp(1)),"_",conds_name(cond_comp(2)),".jpeg");
        saveas(gcf,fullfile(save2,fig_name));
    end
    close;
    close all
end