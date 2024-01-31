function [act_bin] = plot_stat_activity(data,sess_info,stat_cells,phase)
% plot_stat_activity plots activity of stat cells per bin, Z-scored after
% mean over trials
    global num_grp
    if phase == "Task" | phase == "Start"
        clrs = get_color("Enhanced","Suppressed","None");
    else
        clrs = get_color("mix Enhanced","mix Suppressed","mix None");
    end
    % binned
    [act_bin,m_bin_all] = deal(cell(1,length(num_grp)));
    [m_bin,std_bin] = deal(cell(length(num_grp),3));
    [max_m_bin,max_std_bin] = deal(zeros(length(num_grp),3));
%     figure('WindowState','maximized');
    for b = 1:length(num_grp)
        act_bin{b} = cell(1,3);
        for s_i = 1:length(data)
            if ismember(double(sess_info.Session(s_i)),num_grp(b,1):num_grp(b,2))
                for i = 1:size(stat_cells{s_i},2)
                    dat_sess = data{s_i}(logical(stat_cells{s_i}(:,i)),:,:);
                    m_dat_sess = mean(dat_sess,3);
                    zm_dat_sess = zscore(m_dat_sess,[],2);
                    act_bin{b}{i} = [act_bin{b}{i}; zm_dat_sess];
                end
            end
        end
        m_bin_all{b} = cellfun(@(x) max(x(:,25:35),[],2),act_bin{b},'UniformOutput',false);
        m_bin(b,:) = cellfun(@(x) mean(x),act_bin{b},'UniformOutput',false);
        std_bin(b,:) = cellfun(@(x) std(x)/sqrt(size(x,1)),act_bin{b},'UniformOutput',false);
%         subplot(3,2,b)
        figure();
        for ct = 1:size(m_bin,2)
            p(ct) = plot(m_bin{b,ct},'Color',clrs(ct,:),'LineWidth',2);
            hold on
            x2 = [1:length(m_bin{b,ct}), fliplr(1:length(m_bin{b,ct}))];
            inBetween = [(m_bin{b,ct}-std_bin{b,ct}) fliplr((m_bin{b,ct}+std_bin{b,ct}))];
            fill(x2, inBetween,clrs(ct,:),'FaceAlpha',0.3,'EdgeAlpha',0.5);
        end
        ylim([-3 3]);
        yline(30,':k');
%         ylabel('Z-Score');
%         xlabel('Time relative to lick onset(s)');
%         xticks(0:5:45);
%         digits(2);
%         xticklabels(string(vpa(-30:5:15)/30));
%         title(['Bin: ' num2str(b)]);
%         if b == 1
%             legend(p,["Enhanced","Suppressed","None"],'Location','northwest');
%         end
        box off
        % check peak values
        for ct = 1:size(m_bin,2)
            [max_m_bin(b,ct),m_id] = max(m_bin{b,ct});
            max_std_bin(b,ct) = std_bin{b,ct}(m_id);
        end
    end
    figure();
    for ct = 1:size(max_m_bin,2)
        p(ct) = plot(max_m_bin(:,ct),'Color',clrs(ct,:),'LineWidth',2);
        hold on
        x2 = [1:length(max_m_bin(:,ct)), fliplr(1:length(max_m_bin(:,ct)))];
        inBetween = [(max_m_bin(:,ct)-max_std_bin(:,ct))' fliplr((max_m_bin(:,ct)+max_std_bin(:,ct))')];
        fill(x2, inBetween,clrs(ct,:),'FaceAlpha',0.3,'EdgeAlpha',0.5);
    end
%     xlabel("Bin");
%     ylabel("Z-Score");
%     title("Max value per bin, absolute value");
    axis padded
    ylim([0 2.5]);
    xticks(1:length(num_grp));
%     legend(p,["Enhanced","Suppressed","None"],'Location','northwest');
    box off
    % permutation bin1 vs bin6
    num_iters = 5000;
    ct_names = ["Enhanced","Suppressed","None"];
    for i = 1:length(m_bin_all{1})
        pv = permutationTest(m_bin_all{1}{i}, m_bin_all{6}{i},num_iters);
        fprintf('\n %s ,bin 1 vs. bin 6. p-value: %.4f, iters: %d \n',ct_names(i),pv,num_iters);
    end
    [all_dat,cell_id,bin_id] = deal([]);
    for b = 1:length(m_bin_all)
        all_dat = [all_dat ; cell2mat(m_bin_all{b}')];
        cell_id = [cell_id; repelem([1 2 3]',cellfun(@length,m_bin_all{b}),1)];
        bin_id = [bin_id; repelem(b,sum(cellfun(@length,m_bin_all{b})),1)];
    end
    figure();
    [pv,~,stats] = anovan(all_dat,{cell_id,bin_id},'model','interaction','varnames',{'Cell type','Bin'});
    [c,~,~,nms] = multcompare(stats,'Dimension',[1 2],'CType','bonferroni');
    fprintf("post-hoc multiple comparison test bin 1 and 6. p= %.4f",c(15,6))
end