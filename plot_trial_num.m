function plot_trial_num(data,trial_track,sess_info)
    global num_grp
    test_track = cell(1,3);
    sess_track = cell(size(trial_track));
    for trck = 1:length(trial_track)
        if sess_info.Phase(trck) == "Test"
            trial_track{trck}(ismember(trial_track{trck},[1 2]))  = 1;
            trial_track{trck}(ismember(trial_track{trck},[3 4]))  = 2;
            trial_track{trck}(ismember(trial_track{trck},[5 6]))  = 3;
            sess_track{trck} = [sum(trial_track{trck} == 1) sum(trial_track{trck} == 2)];
        end
    end
    for b = 4:length(num_grp)
        test_track{b-3} = [];
        for s_i = 1:length(sess_track)
            if ismember(double(sess_info.Session(s_i)),num_grp(b,1):num_grp(b,2))
                test_track{b-3} = [test_track{b-3} ; (sess_track{s_i}/sum(sess_track{s_i}))*100];
            end
        end
    end
    m_perc = cell2mat(cellfun(@mean ,test_track,'UniformOutput' ,false)');
    mm_perc = mean(m_perc);
    fprintf("Success: Hit: %.2f , Miss: %.2f ",mm_perc(1),mm_perc(2));
    std_perc = cell2mat(cellfun(@(x) std(x)/sqrt(length(x)) ,test_track,'UniformOutput' ,false)');
    clrs = get_color("Hit","Miss");
    figure();
    b = bar(m_perc);
    b(1).FaceColor = clrs(1,:);
    b(2).FaceColor = clrs(2,:);
    hold on
    errorbar([b(1).XEndPoints ;b(2).XEndPoints]',m_perc,std_perc,'k','LineStyle','none');
    xticklabels(4:6);
    legend("Hit","Miss",'Location','northwest');
    box off
    % check how many trials 
    trials_left = cellfun(@(x) size(x,3),data);
    figure();
    subplot(2,2,1); plot(trials_left); title("per session"); 
    xlabel('Session');
    ylabel('Count');
    subplot(2,2,2); histogram(trials_left); title("all");
    xlabel('Trials');
    ylabel('Count');
    sum_conds_test = cellfun(@(x) [sum(ismember(x,[1 2])) sum(ismember(x,[3 4])) ...
        sum(ismember(x,[5 6]))],trial_track(sess_info.Phase == "Test"),'UniformOutput',false);
    sum_conds_train = cellfun(@(x) sum(ismember(x,[1 2])),trial_track(sess_info.Phase == "Train"),'UniformOutput',false);
    subplot(2,2,3); histogram(cellfun(@(x) x(1),sum_conds_train));
    title("training"); 
    xlabel('Trials');
    ylabel('Count');
    X = categorical({'Hit','Miss','FA'});
    X = reordercats(X,{'Hit','Miss','FA'});
    subplot(2,2,4); bar(X,sum(reshape(cell2mat(sum_conds_test),3,[]),2));
    title("test: per condition"); 
    xlabel('Condition');
    ylabel('Count');
    sgtitle('Number of trials per phase');
end