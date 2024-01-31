function [all_data,ind_inv_trials,cnt_licks,licks_time] = LTA_task(TF,train)
% LTA_task crop window around lick in each trial, at task (after sound)
    global wind_pre wind_post End_max_frames twop_fr tosca_fr mouse q_time
    % initialize
    all_data = zeros(size(TF(1).trial_sp,1),wind_pre+wind_post,length(TF));
    ind_inv_trials = zeros(1,length(TF));
    cnt_licks = zeros(1,length(TF));
    time_after_lick = 150; % 5 sec
    licks_time = zeros(length(TF),time_after_lick);
    % go over trials
    for t = 1:length(TF)
        % trial spikes
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
        % check that there is no late licks
        if isempty(lck) | lck(1) + wind_post > size(TF(t).trial_sp,2) 
            all_data(:,:,t) = nan(size(TF(t).trial_sp,1),wind_pre + wind_post);
            ind_inv_trials(t) = 1;
            continue;
        end
        % index of licks after water
        lick_after_w = lck > f_ind(end-1);
        % to count licks
        lck_task = lck(lck >  f_ind(3) & lck <=  f_ind(end));
        if ~isempty(lck_task)
            cnt_licks(t) = 1;
        end
        % lick rate (first lick, train- after water, test - after sound)
        if train
            first_lick_ind = find(lck >= f_ind(end-1));
            if ~isempty(first_lick_ind)
                first_lick_ind = first_lick_ind(1);
                lcs_aligned = lck - lck(first_lick_ind)+1;
                inds = dsearchn(lcs_aligned',[1 time_after_lick]');
                if lcs_aligned(inds(2)) > time_after_lick
                    inds(2) = inds(2) - 1;
                end
                licks_time(t,lcs_aligned(inds(1):inds(2))) = 1;
            end
        else
            first_lick_ind = find(lck >= f_ind(3));
            if ~isempty(first_lick_ind)
                first_lick_ind = first_lick_ind(1);
                lcs_aligned = lck - lck(first_lick_ind)+1;
                inds = dsearchn(lcs_aligned',[1 time_after_lick]');
                if lcs_aligned(inds(2)) > time_after_lick
                    inds(2) = inds(2) - 1;
                end
                licks_time(t,lcs_aligned(inds(1):inds(2))) = 1;
            end
        end
        % take lick that was after quiet X frames and after water
        q_lick = diff(lck) >= q_time; % quiet licks
        after_water_l = lck(lick_after_w);
        % check if first lick is also valid
        first_l = find(lick_after_w);
        if any(first_l)
            first_l = first_l(1);
        end
        if isempty(after_water_l) 
            q_lick = [0 q_lick];
        else
            if first_l == 1 || lck(first_l) -  lck(first_l-1) >= q_time 
                q_lick = [1 q_lick];
            else
                q_lick = [0 q_lick];
            end
        end
        f_lick = lick_after_w == q_lick; % first lick after sound
        f_lick(~lick_after_w) = 0;
        f_lick = lck(f_lick); % frames of valid licks
        if isempty(f_lick) || f_lick(1) + wind_post > size(data,2)
            all_data(:,:,t) = nan(size(TF(t).trial_sp,1),wind_pre + wind_post);
            ind_inv_trials(t) = 1;
            f_lick = [];
        else
            f_lick = f_lick(1); % take only the first valid lick 
            all_data(:,:,t) = data(:,f_lick - wind_pre + 1: f_lick + wind_post);
            f_lick = [];
        end
    end
end