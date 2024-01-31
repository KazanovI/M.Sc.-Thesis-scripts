function [res_fr_bin] = plot_fr_lick(data,sess_info,ind_cells,phase)
    global num_grp
%     FR = cellfun(@(x) mean(squeeze((sum(x,2))/1.5),2),data,'UniformOutput',false); % window is 45 frames which is 1.5 s
%% late addition
    FR1 = cellfun(@(x) mean(squeeze((sum(x(:,1:29,:),2))/0.96),2),data,'UniformOutput',false);
    FR2 = cellfun(@(x) mean(squeeze((sum(x(:,30:45,:),2))/0.53),2),data,'UniformOutput',false);
% binned
    [fr_bin,fr_bin_sess] = deal(cell(2,length(num_grp)));
    [m_bin,std_bin] = deal(zeros(length(num_grp),3,2));
    FRS = [FR1;FR2]
    for t_lick = 1:2
        for b = 1:length(num_grp)
            [fr_bin{t_lick,b},fr_bin_sess{t_lick,b}] = deal(cell(1,3));
            for s_i = 1:length(FRS)
                if ismember(double(sess_info.Session(s_i)),num_grp(b,1):num_grp(b,2))
                    for i = 1:size(ind_cells{s_i},2)
                        fr_bin{t_lick,b}{i} = [fr_bin{t_lick,b}{i}; FRS{t_lick,s_i}(logical(ind_cells{s_i}(:,i)))];
                        fr_bin_sess{t_lick,b}{i} = [fr_bin_sess{t_lick,b}{i} ; mean(FRS{t_lick,s_i}(logical(ind_cells{s_i}(:,i))))];
                    end
                end
            end
            m_bin(b,:,t_lick) = cellfun(@mean ,fr_bin{t_lick,b});
            std_bin(b,:,t_lick) = cellfun(@(x) std(x)/sqrt(length(x)),fr_bin{t_lick,b});
%             res_fr_bin{b} = cell2mat(fr_bin_sess{t_lick,b});
        end
    end
    % errorbar+plot per bin
    if phase == "task" | phase == "start"
        clrs = get_color("Enhanced","Suppressed","None");
    else
        clrs = get_color("mix Enhanced","mix Suppressed","mix None");
    end
    figure();
    subplot(1,2,1)
    for t_lick = 1:2
        for i = 1:size(m_bin,2)
            p(i) = plot(m_bin(:,i,t_lick),'Color',clrs(i,:),'LineWidth',2);
            if t_lick == 2
                p(i).LineStyle = '--';
            end
            hold on
            x2 = [1:length(m_bin(:,i,t_lick)), fliplr(1:length(m_bin(:,i,t_lick)))];
            inBetween = [(m_bin(:,i,t_lick)-std_bin(:,i,t_lick))' fliplr((m_bin(:,i,t_lick)+std_bin(:,i,t_lick))')];
            fill(x2, inBetween,clrs(i,:),'FaceAlpha',0.3,'EdgeAlpha',0.5);
        end
        title("Mean FR. pre: - , lick: --");
        xlabel("Bin");
        ylabel("Mean FR(Hz)");
        axis padded
    %     legend(p,["Enhanced","Suppressed","None"]);
        box off
    end
    subplot(1,2,2)
    m_bin_ratio = m_bin(:,:,2)./m_bin(:,:,1);
    std_bin_ratio = std_bin(:,:,2)./std_bin(:,:,1);
    for i = 1:size(m_bin_ratio,2)
        p(i) = plot(m_bin_ratio(:,i),'Color',clrs(i,:),'LineWidth',2);
        hold on
%         x2 = [1:length(m_bin_ratio(:,i)), fliplr(1:length(m_bin_ratio(:,i)))];
%         inBetween = [(m_bin_ratio(:,i)-std_bin_ratio(:,i))' fliplr((m_bin_ratio(:,i)+std_bin_ratio(:,i))')];
%         fill(x2, inBetween,clrs(i,:),'FaceAlpha',0.3,'EdgeAlpha',0.5);
    end
    title("mean FR ratio");
    xlabel("Bin");
    ylabel("Ratio");
    axis padded
    
    
    %%
    % binned
    [fr_bin,fr_bin_sess] = deal(cell(1,length(num_grp)));
    [m_bin,std_bin] = deal(zeros(length(num_grp),3));
%     cells_bin = cell(1,length(num_grp));
    for b = 1:length(num_grp)
        [fr_bin{b},fr_bin_sess{b}] = deal(cell(1,3));
%         cells_bin{b} = 0;
        for s_i = 1:length(FR)
            if ismember(double(sess_info.Session(s_i)),num_grp(b,1):num_grp(b,2))
                for i = 1:size(ind_cells{s_i},2)
                    fr_bin{b}{i} = [fr_bin{b}{i}; FR{s_i}(logical(ind_cells{s_i}(:,i)))];
                    fr_bin_sess{b}{i} = [fr_bin_sess{b}{i} ; mean(FR{s_i}(logical(ind_cells{s_i}(:,i))))];
                end
%                 cells_bin{b} = cells_bin{b} + size(data{b},1);
            end
        end
        m_bin(b,:) = cellfun(@mean ,fr_bin{b});
        std_bin(b,:) = cellfun(@(x) std(x)/sqrt(length(x)),fr_bin{b});
        res_fr_bin{b} = cell2mat(fr_bin_sess{b});
    end
    % errorbar+plot per bin
    if phase == "task" | phase == "start"
        clrs = get_color("Enhanced","Suppressed","None");
    else
        clrs = get_color("mix Enhanced","mix Suppressed","mix None");
    end
    figure();
    for i = 1:size(m_bin,2)
        p(i) = plot(m_bin(:,i),'Color',clrs(i,:),'LineWidth',2);
        hold on
        x2 = [1:length(m_bin(:,i)), fliplr(1:length(m_bin(:,i)))];
        inBetween = [(m_bin(:,i)-std_bin(:,i))' fliplr((m_bin(:,i)+std_bin(:,i))')];
        fill(x2, inBetween,clrs(i,:),'FaceAlpha',0.3,'EdgeAlpha',0.5);
    end
%     title("Mean FR per bin and cell type, lick time");
%     xlabel("Bin");
%     ylabel("Mean FR(Hz)");
    axis padded
%     ylim([20 60]);
%     legend(p,["Enhanced","Suppressed","None"]);
    box off
    % statistics ANOVA2
    [all_dat,cell_id,bin_id] = deal([]);
    for b = 1:length(fr_bin)
        all_dat = [all_dat ; cell2mat(fr_bin{b}')];
        cell_id = [cell_id; repelem([1 2 3]',cellfun(@length,fr_bin{b}),1)];
        bin_id = [bin_id; repelem(b,sum(cellfun(@length,fr_bin{b})),1)];
    end

    [pv,~,stats] = anovan(all_dat,{cell_id,bin_id},'model','interaction','varnames',{'Cell type','Bin'});
    [c,~,~,nms] = multcompare(stats,'Dimension',[1 2]);
    fprintf("post-hoc multiple comparison test bin 1 and 4. p= %.4f",c(9,6))
    % permutations bin 1 vs. 4
    num_iters = 5000;
    ct_names = ["Enhanced","Suppressed","None"];
    for i = 1:length(fr_bin{1})
        pv = permutationTest(fr_bin{1}{i}, fr_bin{4}{i},num_iters);
        fprintf('\n %s ,bin 1 vs. bin 4. p-value: %.4f, iters: %d \n',ct_names(i),pv,num_iters);
    end
    % histograms
    figure();
    for b = 1:length(fr_bin)
        ax(b) = subplot(2,3,b);
        for clt = 1:length(fr_bin{b})
            xtk = [min(fr_bin{b}{clt}):5:max(fr_bin{b}{clt})];
            perc = (histcounts(fr_bin{b}{clt},length(xtk))/length(fr_bin{b}{clt}))*100;
            br(clt) = bar(perc,'FaceColor',clrs(clt,:),'FaceAlpha',0.5);
            bmean(clt) = mean(fr_bin{b}{clt});
            hold on
        end
        title(['Bin: ' num2str(b), 'Mean: ' num2str(bmean)]);
        xlabel('FR(Hz)');
        ylabel('% of cells');
        xticks(1:2:length(xtk));
        xticklabels(xtk(1:2:end));
    end
    legend(br,["Enhanced","Suppressed","None"]);
    linkaxes(ax,'y');
    
    % by mouse
    mice = unique(sess_info.Mouse);
    figure();
    for m = 1:length(mice)
        fr_bin = cell(1,length(num_grp));
        [m_bin,std_bin] = deal(zeros(length(num_grp),3));
        for b = 1:length(num_grp)
            fr_bin{b} = cell(1,3);
            for s_i = 1:length(sess_info.Mouse)
                if sess_info.Mouse(s_i) == mice(m) & ismember(double(sess_info.Session(s_i)),num_grp(b,1):num_grp(b,2))
                    for i = 1:size(ind_cells{s_i},2)
                        fr_bin{b}{i} = [fr_bin{b}{i}; FR{s_i}(logical(ind_cells{s_i}(:,i)))];
                    end
                end
            end
            m_bin(b,:) = cellfun(@mean ,fr_bin{b});
            std_bin(b,:) = cellfun(@(x) std(x)/sqrt(length(x)),fr_bin{b});
        end
        subplot(3,3,m)
        for i = 1:size(m_bin,2)
            p(i) = plot(m_bin(:,i),'Color',clrs(i,:),'LineWidth',2);
            hold on
            x2 = [1:length(m_bin(:,i)), fliplr(1:length(m_bin(:,i)))];
            inBetween = [(m_bin(:,i)-std_bin(:,i))' fliplr((m_bin(:,i)+std_bin(:,i))')];
            inBetween(isnan(inBetween)) = 0;
            fill(x2, inBetween,clrs(i,:),'FaceAlpha',0.3,'EdgeAlpha',0.5);
        end
        title(mice(m));
        xlabel("Bin");
        ylabel("Mean FR(Hz)");
        axis padded
        if m == length(mice)
            legend(p,["Enhanced","Suppressed","None"]);
        end
        box off
    end
end