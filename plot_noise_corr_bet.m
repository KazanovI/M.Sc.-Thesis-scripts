function plot_noise_corr_bet(data,sess_info,inds,phase)
    global num_grp
    window = [15 29;30 45]; % before and after lick
    mrk = ['d','x'];
    tit = ["Enh-Supp","Enh-None","Supp-None"];
    cmps = [1,2;1,3;2,3];
    [num_pairs,dat_pairs] = deal(cell(1,length(tit)));
    figure();
    for ct = 1:length(tit)
        for wind = 1:length(window)
            [mean_bin,std_bin] = deal(cell(1,length(num_grp)));
            for b = 1:size(num_grp,1) % per bin
                bin_cor = []; % initialize
                for s_i = 1:length(data)% session
                    % check if in bin
                    if ismember(double(sess_info.Session(s_i)),num_grp(b,1):num_grp(b,2))
                        dat_cell = [{data{s_i}(logical(inds{s_i}(:,cmps(ct,1))),:,:)}...
                            {data{s_i}(logical(inds{s_i}(:,cmps(ct,2))),:,:)}];
                        % check both cell types are valid
                        if any(cellfun(@(x) all(isnan(x),'all'),dat_cell))
                            continue;
                        else
                            if length(window(wind,1):window(wind,2)) ~= size(dat_cell{1},2)
                                dat_cell = cellfun(@(x) x(:,window(wind,1):window(wind,2),:),dat_cell,'UniformOutput',false);
                            end
                            % calc pair-wise correlations, first output is
                            % mean over trials,second is all corrs
                            [bet_cor,bet_cor_all] = calc_bet_corr2(dat_cell{1},dat_cell{2});
                        end
                        bin_cor = [bin_cor ; bet_cor];
                    end
                end
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
                errorbar(b+locs(wind),mean_bin{b},std_bin{b},'Color','k','Marker',mrk(wind));
            end
        end
        axis('padded')
        xlim([0 length(num_grp)+1])
        xticks(1:length(num_grp));
        xticklabels(1:6);
        xlabel("Bin");
        ylabel("r");
        title(['Noise corr, cells: ' char(tit(ct))]);
    end
    linkaxes(ax,'y');
    sgtitle([char(phase) ', between groups: \diamondsuit: before, \times: Lick']);

    % plot differences (after lick - before)
    diff_pairs = cell(1,length(dat_pairs));
    figure();
    for ct = 1:length(dat_pairs)
        ax(ct) = subplot(3,1,ct);
        for b = 1:length(dat_pairs{ct}{1})
            diff_pairs{ct}{b} = dat_pairs{ct}{2}{b} - dat_pairs{ct}{1}{b};
            hold on;
            errorbar(b,mean(diff_pairs{ct}{b},'omitnan'),std(diff_pairs{ct}{b},[],1,'omitnan')/sqrt(size(diff_pairs{ct}{b},1))...
                ,'Color','k','Marker','o');
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
    sgtitle([char(phase) ', Between groups: difference of noise corr, after lick - before']);

    % plot train and test separate
    if phase == "task" | phase == "start"
        clrs = get_color("Enhanced","Suppressed","None");
        clrs = clrs * 0.9;
    else
        clrs = get_color("mix Enhanced","mix Suppressed","mix None");
        clrs = clrs * 0.9;
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
        xticks(1:3);
        xticklabels(bins_id(b,:));
        xlabel('Bin');
        ylabel('Diff corr');
        legend(pp,tit,'Location','northeastoutside');
        title([char(phase) ', Between groups, ' char(tit_time(b))]);
    end
    linkaxes(ax,'y');
end