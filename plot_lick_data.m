function [ret_all_res] = plot_lick_data(licks_count,lcks_rate,valid_sess_LTA,trial_track)
% count - both start and task data
% rate - only task, lick frames
    global num_grp
    eng = zeros(1,length(licks_count));
    cnt_both = zeros(2,length(licks_count));
    for ss = 1:length(licks_count)
        trials_ind = round(prctile(1:length(licks_count{1,ss}),[0 100]));
        trials_ind = trials_ind(1):trials_ind(2);
        cnt_both(:,ss) = cellfun(@(x) (sum(x(trials_ind))/length(trials_ind))*100,licks_count(:,ss));
        % calc engagment, lick atleast in start or task
        eng(ss) = (1-sum(sum(cell2mat(licks_count(:,ss))) == 0)/length(licks_count{1,ss}));
    end
    % binned
    [all_res,eng_bin] = deal(cell(1,length(num_grp)));
    for b = 1:length(num_grp)
        for si = 1:length(cnt_both)
            if ismember(str2num(valid_sess_LTA.Session(si)),[num_grp(b,1):num_grp(b,2)])
                all_res{1,b} = [all_res{1,b} , cnt_both(:,si)];
                eng_bin{b} = [eng_bin{b} , eng(si)];
            end
        end
    end
    ret_all_res = all_res;
    figure();
    for b = 1:length(all_res)
    %     s = scatter(repelem([b-0.2 b+0.2],size(all_res{b},2),1),all_res{b}',[],[1 0 0; 0 0 1]);
    %     hold on;
        s = bar(b,mean(all_res{b}'),'BarWidth',1);
        s(1).FaceColor = [0 0 0];
        s(2).FaceColor = [1 1 1];
        endpnts = [s.XEndPoints];
        hold on;
        errorbar(endpnts,mean(all_res{b}'),std(all_res{b}')/sqrt(size(all_res{b},2)),'k+');
    end
    xlim([0 7]);
    ylim([0 100]);
    xticks(1:6);
    xticklabels(1:6);
    xlabel('Bin');
    ylabel('% licks ');
    legend(s,["Start","Task"]);
    title("% of trials with lick by start or task");
    box off
    % 2way anova (start/task & train/test)
    % construct matrix
    [dat_vec,time,state,cmp_bin,dat_bin] = deal([]); % data vectorized,start/task,train/test
    for i = 1:length(all_res)
        dat_temp = all_res{i}';
        time = [time; repelem(["start","task"],size(dat_temp,1))'];
        cmp_bin = [cmp_bin; repelem(i,size(dat_temp,1))'];
        dat_bin = [dat_bin; all_res{i}(2,:)'];
        if i<=3
            state = [state; repelem(["train"],size(dat_temp,1)*2)'];
        else
            state = [state; repelem(["test"],size(dat_temp,1)*2)'];
        end
        dat_vec = [dat_vec; dat_temp(:)];
    end
    p = anovan(dat_vec,{time state},'model',1,'varnames',{'time','state'});
    [~,~,stats] = anova1(dat_bin,cmp_bin); % anova for lick in task between bins
    multcompare(stats);
    % engagment
    figure();
    for b = 1:length(all_res)
        s = scatter(repelem(b,size(eng_bin{b},2),1),eng_bin{b},25,[0.3 0.3 0.3],'filled'...
            ,'jitter','on','jitterAmount',0.2,'MarkerFaceAlpha',0.8);
        hold on;
        errorbar(b,mean(eng_bin{b}),std(eng_bin{b})/sqrt(size(eng_bin{b},2)),'ko',...
            'MarkerFaceColor','auto','MarkerSize',5);
    end
    xlim([0 7]);
    xticks(1:6);
    ylim([0.2 1])
    xlabel('Bin');
    ylabel('Engagment');
    title('% of trials with lick(any state)');
    % 1way anova (bin)
    % construct matrix
    [cmp_bin,dat_bin] = deal([]); % data vectorized,start/task,train/test
    for i = 1:length(eng_bin)
        cmp_bin = [cmp_bin; repelem(i,size(dat_temp,1))'];
        dat_bin = [dat_bin; all_res{i}(2,:)'];
    end
    [~,~,stats] = anova1(dat_bin,cmp_bin);
    multcompare(stats);

    % plot lick rate by mouse and session
    mice = unique(valid_sess_LTA.Mouse);
    for m = 1:length(mice)
        figure();
        ind_sess = find(valid_sess_LTA.Mouse == mice(m));
        nrows = floor(sqrt(length(ind_sess)));
        ncols = ceil(length(ind_sess)/nrows);
        for si = 1:length(ind_sess)
            subplot(nrows, ncols, si);
            plot(mean(lcks_rate{ind_sess(si)}));
            title(append("Session: ", valid_sess_LTA.Session(ind_sess(si))));
        end
        sgtitle(valid_sess_LTA.Mouse(ind_sess(1)));
    end
    %% also per condition
    mice = unique(valid_sess_LTA.Mouse);
    cnds = [1 2;3 4;5 6];
    clrs = get_color("Hit","Miss","FA");
    lcks_by_cond = cell(1,length(mice));
    lcks_train = cell(1,length(mice));
    for m = 1:length(mice)
        figure();
        ind_sess = find(valid_sess_LTA.Mouse == mice(m));
        nrows = floor(sqrt(length(ind_sess)));
        ncols = ceil(length(ind_sess)/nrows);
        for si = 1:length(ind_sess)
            subplot(nrows, ncols, si);
            sess_cond = trial_track{ind_sess(si)};
            if length(unique(sess_cond)) > 2
                for cn = 1:size(cnds,1)
                    lcks_by_cond{m}{cn,si} = mean(lcks_rate{ind_sess(si)}(ismember(sess_cond,cnds(cn,:)),:));
                    if cn == 1
                        lck_hit = lcks_rate{ind_sess(si)}(ismember(sess_cond,cnds(cn,:)),:);
%                         lck_hit_shift = circshift(lck_hit,-30,1);
                        plot(mean(lck_hit),'Color',clrs(cn,:),'LineWidth',2);
                    else
                        plot(mean(lcks_rate{ind_sess(si)}(ismember(sess_cond,cnds(cn,:)),:)),'Color',clrs(cn,:),'LineWidth',2);
                    end
                    
                    hold on;
                end
            else
                lcks_train{m}{si} = mean(lcks_rate{ind_sess(si)});
                plot(mean(lcks_rate{ind_sess(si)}));
            end
            title(append("Session: ", valid_sess_LTA.Session(ind_sess(si))));
        end
        sgtitle(valid_sess_LTA.Mouse(ind_sess(1)));
    end
    lcks_by_cond = lcks_by_cond(cellfun(@(x) ~isempty(x),lcks_by_cond));
    per_mice = cell(1,3);
    for m = 1:length(lcks_by_cond)
        dat_m = lcks_by_cond{m}(:,cellfun(@(x) ~isempty(x),lcks_by_cond{m}(1,:)));
        for cnd = 1:size(dat_m,1)
            per_mice{cnd} = [per_mice{cnd} ; cell2mat(dat_m(cnd,:)')];
        end
    end
    per_mice_s = cellfun(@(x) smoothdata(x,2,'movmean',[0 3]),per_mice,'UniformOutput',false);
    for cnd = 1:length(per_mice_s)
        per_mice_s{cnd}(:,1) = per_mice{cnd}(:,1);
    end
    m_all = cellfun(@(x) mean(x),per_mice_s,'UniformOutput',false);
    sd_all = cellfun(@(x) std(x)/sqrt(size(x,1)),per_mice_s,'UniformOutput',false);
    figure();
    for i = 1:length(m_all)
        p(i) = plot(m_all{i},'Color',clrs(i,:),'LineWidth',2);
        hold on 
        x2 = [1:length(m_all{i}), fliplr(1:length(m_all{i}))];
        inBetween = [(m_all{i}-sd_all{i}) fliplr((m_all{i}+sd_all{i}))];
        fill(x2, inBetween,clrs(i,:),'FaceAlpha',0.3,'EdgeAlpha',0.5);
    end
    xticks(0:20:length(m_all{1}));
%     digits(2);
%     xticklabels(string(vpa(0:20:length(m_all{1}))/30));
    set(gca,'XScale','log')
    xlabel("Time relative to lick onset(s)");
    ylabel("Proportion of trials");
    title("Licks by condition");
    legend(p,["Hit","Miss","FA"]);
    box off
end