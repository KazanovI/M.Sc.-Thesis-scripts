function plot_enh_peak(task,mix)
    global num_grp
    nms = ["Enhanced","Suppressed","None"];
    both = [task;mix];
    for ct = 1:length(nms)
        for b = 1:6
            dat = cellfun(@(x) x{ct},both(:,b),'UniformOutput',false);
            if ct == 2
                dat = cellfun(@(x) -x,dat,'UniformOutput',false);
            end
%             m_bin_all{ct}(:,b) = cellfun(@(x) max(x(:,25:35),[],2),dat,'UniformOutput',false);
            m_bin(:,b) = cellfun(@(x) mean(x),dat,'UniformOutput',false);
            std_bin(:,b) = cellfun(@(x) std(x)/sqrt(size(x,1)),dat,'UniformOutput',false);
            % check peak values
            [max_m_bin(:,b),m_id] = cellfun(@(x) max(x),m_bin(:,b));
            for ti = 1:2
                max_std_bin(ti,b) = std_bin{ti,b}(m_id(ti));
                max_dat_all{ct}{b}(:,ti) = dat{ti}(:,m_id(ti));
            end
        end
        nms_ct = [nms(ct) append("mix ",nms(ct))];
        clrs = get_color(nms_ct);
        figure();
        for i = 1:size(max_m_bin,1)
            p(i) = plot(max_m_bin(i,:),'Color',clrs(i,:),'LineWidth',2);
            hold on
            x2 = [1:length(max_m_bin(i,:)), fliplr(1:length(max_m_bin(i,:)))];
            inBetween = [(max_m_bin(i,:)-max_std_bin(i,:))'; flipud((max_m_bin(i,:)+max_std_bin(i,:))')];
            fill(x2, inBetween,clrs(i,:),'FaceAlpha',0.3,'EdgeAlpha',0.5);
        end
%         xlabel("Bin");
%         ylabel("Z-Score");
    %     title("Max value per bin, absolute value");
        axis padded
        ylim([0 2.5]);
        xticks(1:length(num_grp));
%         legend(p,nms_ct,'Location','northwest');
        box off
        % stat
%         num_iters = 5000;
%         % bin 6 comparison
%         pv6 = permutationTest(max_dat_all{6}(:,1), max_dat_all{6}(:,2),num_iters);
%         % bin 1vs6 comparison
%         pv16task = permutationTest(max_dat_all{1}(:,1), max_dat_all{6}(:,1),num_iters);
%         pv16mix = permutationTest(max_dat_all{1}(:,2), max_dat_all{6}(:,2),num_iters);
        [all_dat,time_id,bin_id] = deal([]);
        for b = 1:length(max_dat_all{ct})
            all_dat = [all_dat ;  max_dat_all{ct}{b}(:)];
            time_id = [time_id; repelem([1 2]',length(max_dat_all{ct}{b}))];
            bin_id = [bin_id; repelem(b,length(max_dat_all{ct}{b})*2)'];
        end
        figure();
        [pv,tbl,stats] = anovan(all_dat,{time_id,bin_id},'model','interaction','varnames',{'Time','Bin'});
        [c,~,~,nms2] = multcompare(stats,'Dimension',[1 2],'CType','bonferroni');
        fprintf("Cell: %d bin 6 between times. p= %.4f \n",ct,c(66,6));
        fprintf("Cell: %d bin 1vs6 task|mix. p= %.4f|%.4f \n",ct,c(10,6),c(21,6));
    end
    % anova3
%     [all_dat,cell_id,bin_id,time_id] = deal([]);
%     for ct = 1:length(m_bin_all)
%         for b = 1:length(m_bin_all{ct})
%             all_dat = [all_dat ; cat(1,m_bin_all{ct}{:,b})];
%             cell_id = [cell_id; repelem(ct,length(cat(1,m_bin_all{ct}{:,b})))'];
%             bin_id = [bin_id; repelem(b,length(cat(1,m_bin_all{ct}{:,b})))'];
%             time_id = [time_id; repelem([1 2]',length(m_bin_all{ct}{1,b}))]; 
%         end
%     end
%     figure();
%     [pv,~,stats] = anovan(all_dat,{cell_id,bin_id,time_id},'model','interaction','varnames',{'Cell type','Bin','Time'});
%     [c,~,~,nms] = multcompare(stats,'Dimension',[1 2 3],'CType','bonferroni');
% close all
end