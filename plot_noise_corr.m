function plot_noise_corr(data,sess_info,inds,trial_track,phase)
    global num_grp
    window = [1:15;25:39]; % before lick
    % within stat cell group
    tit = ["Enhanced","Suppressed","None"];
    wnd_clr = ["w","k"];
    wnd_name = ["Before","During"];
    [num_pairs,dat_pairs,dat_pairs_sess] = deal(cell(1,length(tit)));
    figure();
    for ct = 1:length(tit)
        for wind = 1:size(window,1)
            [mean_bin,std_bin] = deal(cell(1,length(num_grp)));
            for b = 1:size(num_grp,1) % per bin
                bin_cor = []; % initialize
                for s_i = 1:length(data)% session
                    % check if in bin
                    if ismember(double(sess_info.Session(s_i)),num_grp(b,1):num_grp(b,2))
                        dat_cell = data{s_i}(logical(inds{s_i}(:,ct)),:,:);
                        if size(dat_cell,1) < 3 || all(isnan(dat_cell),'all') 
                            continue;
                        else
                            if length(window) ~= size(dat_cell,2)
                                dat_cell = dat_cell(:,window(wind,:),:);
                            end
                            % calc pair-wise correlations, first output is
                            % mean over trials,second is all corrs
                            [bet_cor,~,~] = calc_bet_corr2(dat_cell);
                        end
                        bin_cor = [bin_cor ; bet_cor];
                        dat_pairs_sess{ct}{b}{wind}{s_i} = bet_cor;
                    end
                end
                % limit max pair number
                num_pairs{ct}{b} = sum(~isnan(bin_cor));
                max_pairs = inf;
                if size(bin_cor,1) >  max_pairs
                    ind_cells = randperm(size(bin_cor,1),max_pairs);
                else
                    ind_cells = 1:size(bin_cor,1);
                end
                dat_pairs{ct}{wind}{b} = bin_cor(ind_cells,:);
                mean_bin{b} = mean(bin_cor(ind_cells,:),1,'omitnan'); % mean over cells
                std_bin{b} = std(bin_cor(ind_cells,:),[],1,'omitnan')/sqrt(size(bin_cor(ind_cells,:),1));
            end
            ax(ct) = subplot(3,1,ct);
            for b = 1:length(num_grp)
                hold on;
                locs = [-0.2,0.2];
                errorbar(b+locs(wind),mean_bin{b},std_bin{b},'Color','k');
                if b < length(num_grp)
                    bar(b+locs(wind),mean_bin{b},'FaceColor',wnd_clr(wind),'FaceAlpha',0.5,...
                        'BarWidth',0.3);
                else
                    br(wind) = bar(b+locs(wind),mean_bin{b},'FaceColor',wnd_clr(wind),'FaceAlpha',0.5,...
                        'BarWidth',0.3,'DisplayName',wnd_name(wind));
                end
            end
        end
        legend(br);
        axis('padded')
        xlim([0 length(num_grp)+1])
        xticks(1:length(num_grp));
        xticklabels(1:6);
        xlabel("Bin");
        ylabel("r");
        title(['Noise corr, cells: ' char(tit(ct))]);
    end
    linkaxes(ax,'y');
    sgtitle('In group');
    
    % plot corrs per bin and cell type and time
    clrs2 = get_color(["Enhanced","Suppressed","None"]);
%     figure();
    for ct = 1:length(dat_pairs)
%         ax(ct) = subplot(length(dat_pairs),1,ct);
        figure();
        for tim = 1:length(dat_pairs{ct})
            hold on;
            mns = cellfun(@(x) mean(x,'omitnan'),dat_pairs{ct}{tim});
            stds = cellfun(@(x) std(x,'omitnan')/sqrt(length(x)) ,dat_pairs{ct}{tim});
            p(tim) = plot(mns,'Color',clrs2(ct,:),'LineWidth',2);
            if tim == 1
                p(tim).LineStyle = '-';
            else
                p(tim).LineStyle = ':';
            end
            x2 = [1:length(mns), fliplr(1:length(mns))];
            inBetween = [(mns-stds) fliplr(mns+stds)]; % vector,same length as x2
            fill(x2, inBetween,clrs2(ct,:),'FaceAlpha',0.3,'EdgeAlpha',0.5);
        end
        xlim([0 length(num_grp)+1])
        xticks(1:length(num_grp));
        xticklabels(1:6);
        xlabel("Bin");
        ylabel("r");
        ylim([0 0.1]);
%         title(['Noise corr, cells: ' char(tit(ct))]);
%         if ct == 1
            legend(p,"Before","During");
%         end
        %ANOVA2
        time = [];
        bin = [];
        dat_all = [];
        len_bins = cellfun(@length, dat_pairs{ct}{1});
        for tt = 1:length(dat_pairs{ct})
            time = [time ; repelem(tt,sum(len_bins),1)];
            bin = [bin ; repelem([1:6]',len_bins,1)];
            dat_all = [dat_all ; cell2mat(dat_pairs{ct}{tt}')];
        end
        figure();
        [~,tbl,stats] = anovan(dat_all,{time bin},'model',2,'varnames',{'time','Bin'});
        [c,~,~,gnms] = multcompare(stats,'Dimension',[1 2],'CType','bonferroni');
        ns = string(gnms(c(c(:,6) > 0.05,1:2)));
        for n = 1:size(ns,1)
            fprintf('NS %s: %s and %s \n',string(ct),ns(n,1),ns(n,2));
        end
    end
%     linkaxes(ax,'y');
    



    % plot differences (after lick - before)
    diff_pairs = cell(1,length(dat_pairs));
    figure();
    for ct = 1:length(dat_pairs)
        ax(ct) = subplot(3,1,ct);
        for b = 1:length(dat_pairs{ct}{1})
            diff_pairs{ct}{b} = dat_pairs{ct}{2}{b} - dat_pairs{ct}{1}{b};
            hold on;
            errorbar(b,mean(diff_pairs{ct}{b},'omitnan'),std(diff_pairs{ct}{b},[],1,'omitnan')/sqrt(size(diff_pairs{ct}{b},1))...
                ,'Color','k','Marker','x');
        end
        axis('padded');
        xlim([0 length(num_grp)+1]);
        xticks(1:length(num_grp));
        xticklabels(1:6);
        xlabel("Bin");
        ylabel("diff r");
        title(['Noise corr, cells: ' char(tit(ct))]);
    end
    linkaxes(ax,'y');
    sgtitle('In group, difference of noise corr, after lick - before');
    
    % plot train and test separate
    if phase == "task" | phase == "start"
        clrs = get_color("Enhanced","Suppressed","None");
    else
        clrs = get_color("mix Enhanced","mix Suppressed","mix None");
    end
    bins_id = [1:3;4:6];
    tit_time = ["Train","Test"];
    for b = 1:size(bins_id,1)
        figure();
        ax(b) = axes;
        hold on
        for ct = 1:length(diff_pairs)
            m_r = cellfun(@(x) mean(x,'omitnan'),diff_pairs{ct}(bins_id(b,:)));
            m_std = cellfun(@(x) std(x,'omitnan')/sqrt(length(x)),diff_pairs{ct}(bins_id(b,:)));
            pp(ct) = plot(m_r,'LineWidth',2,'Color',clrs(ct,:));
            x2 = [1:length(m_r), fliplr(1:length(m_r))];
            inBetween = [(m_r-m_std) fliplr(m_r+m_std)]; % vector,same length as x2
            fill(x2, inBetween,clrs(ct,:),'FaceAlpha',0.3,'EdgeAlpha',0.5);
        end
        axis padded
        ylim([-0.015 0.05]);
        xticks(1:3);
        xticklabels(bins_id(b,:));
        xlabel('Bin');
        ylabel('Diff corr');
        legend(pp,tit,'Location','best');
        title(['In group, ' char(tit_time(b))]);
        box off
    end
%     linkaxes(ax,'y');
    
    % distribution of corrs train-test
    bins_id = [1 3;4 6];
    for b = 1:size(bins_id,1)
        figure();
        for ct = 1:length(diff_pairs)
            ax(ct) = subplot(length(diff_pairs),1,ct);
            histogram(diff_pairs{ct}{bins_id(b,1)},'Normalization','pdf','BinMethod','scott');
            hold on
            histogram(diff_pairs{ct}{bins_id(b,2)},'Normalization','pdf','BinMethod','scott');
            [~,p] = kstest2(diff_pairs{ct}{bins_id(b,1)},diff_pairs{ct}{bins_id(b,2)});
            title({[char(tit(ct))],['p: ' num2str(p)]});
            xlabel('Diff Corr');
            ylabel('Probability');
            if ct == 1
                legend(num2str(bins_id(b,1)),num2str(bins_id(b,2)));
            end
            box off
        end
        sgtitle({['PDF of corrs, In group, ' char(tit_time(b))] ,['Sessions: ' num2str(bins_id(b,1)) ' VS. ' num2str(bins_id(b,2))]});
        linkaxes(ax,'xy');
    end
    
    % distribution of corrs over bins
%     disb = cell(1,length(dat_pairs)); % to compare dists
%     xclr = ['b','r'];
%     figure();
%     for ct = 1:length(dat_pairs)
%         ax(ct) = subplot(length(dat_pairs),1,ct);
%         for wind = 1:length(dat_pairs{ct})
%             h_dat = cell2mat(dat_pairs{ct}{wind}');
%             [count,edg] = histcounts(h_dat,[-0.3:0.005:0.3]);
%             disb{ct}{wind} = (count/length(h_dat))*100;
%             br(wind) = bar((count/length(h_dat))*100);
%             br(wind).FaceAlpha = 0.5;
%             [~,loc] = min(abs(edg - mean(h_dat,'omitnan')));
%             xline(loc,xclr(wind));
%             hold on;
%         end
%         title(char(tit(ct)));
%         ylabel('percentage(%)');
%         xlabel('corr');
%         xticks(1:5:length(edg));
%         xticklabels(edg(1:5:end));
%         if ct==1
%             legend(br,"Before","After");
%         end
%         box off
%     end
%     sgtitle('In group, distribution of differences per cell type');
%     linkaxes(ax,'y');
end