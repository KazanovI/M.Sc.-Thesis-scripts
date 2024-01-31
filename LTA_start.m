function [all_data,ind_inv_trials,cnt_licks] = LTA_start(TF)
% LTA_start crop window around lick in each trial, at start (before sound)
    global wind_pre wind_post End_max_frames twop_fr tosca_fr mouse q_time
    all_data = zeros(size(TF(1).trial_sp,1),wind_pre+wind_post,length(TF));
    all_data(:,:,1) = nan;
    ind_inv_trials = zeros(1,length(TF));
    cnt_licks = zeros(1,length(TF));
    % go over trials
    for t = 1:length(TF)
        if size(TF(t).States.End.sp,2) < End_max_frames
            data = TF(t).trial_sp;
        else
            data = TF(t).trial_sp(:,1:End_max_frames);
        end
        % get fields data
        f_names = string(fieldnames(TF(t).States));
        for fn = 1:length(f_names)
            f_ind(fn) = ceil(TF(t).States.(f_names{fn}).Tosca_Frames(1)*(twop_fr/tosca_fr));
        end
        % get licks
        if any(TF(t).Lick)
            lck = unique(ceil(find(TF(t).Lick)*(twop_fr/tosca_fr)));
        elseif any(contains(string(fieldnames(TF(t))),"lickometer_2"))
            lck = unique(ceil(find(TF(t).lickometer_2)*(twop_fr/tosca_fr)));
        else
            lck = unique(ceil(find(TF(t).lickometer_1)*(twop_fr/tosca_fr)));
        end
        % to count licks
        lck_start = lck(lck <=  f_ind(3));
        if ~isempty(lck_start)
            cnt_licks(t) = 1;
        end
        % check that there is no late licks
        if isempty(lck) | lck(1) + wind_post > size(TF(t).trial_sp,2) 
            all_data(:,:,t) = nan(size(TF(t).trial_sp,1),wind_pre + wind_post);
            ind_inv_trials(t) = 1;
            continue;
        end
        % section sta start
        lick_bef_w = lck < f_ind(2);
        if t>1 & any(lick_bef_w)
            before_water_l = lck(lick_bef_w); % inidices of licks in current trial
            q_lick = [0 diff(before_water_l) >= q_time];% quiet lick in currnet trial
            data_prev = TF(t-1).trial_sp; 
            lck_bef_trial = unique(ceil(find(TF(t-1).Lick(end-(tosca_fr-1):end))*(twop_fr/tosca_fr))); % prev trial licks
            if before_water_l(1) >= q_time
                q_lick(1) = 1;
            end
            f_lick_bef = before_water_l(logical(q_lick));
            data_both = [data_prev data];
            new_ind = size(data_prev,2)+f_lick_bef;
            % if lick in current trial is not "quiet", check prev trial licks
            if any(q_lick) || any(before_water_l(1)==1:q_time-1) & ~any(lck_bef_trial > twop_fr - abs(q_time - before_water_l(1)))
                if isempty(f_lick_bef)
                    f_lick_bef = before_water_l(1);
                    new_ind = size(data_prev,2)+f_lick_bef;
                end
                f_lick_bef = f_lick_bef(1);
                % save lick data
                all_data(:,:,t) = data_both(:,new_ind - wind_pre + 1: new_ind + wind_post);
            else
                % nan
                all_data(:,:,t) = nan(size(TF(t).trial_sp,1),wind_pre + wind_post);
                f_lick_bef = [];
                ind_inv_trials(t) = 1;
            end
        else
            % nan
            all_data(:,:,t) = nan(size(TF(t).trial_sp,1),wind_pre + wind_post);
            f_lick_bef = [];
            ind_inv_trials(t) = 1;
        end
    end
end