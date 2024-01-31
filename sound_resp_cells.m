function [ret_responsive_cells] = sound_resp_cells(TF)
% all trials - sound responsiveness
%parameters
global tosca_fr twop_fr
    roi_wind = 10; % where to look peaks at
    sound_wind = 1; % frames of window of baseline & sound response 
    method = 2; % method to choose peak frames
    split = 1; % divide by trials with same sound
    by_sound = zeros(2,length(TF));
    tf_start = cell(1,length(TF));
    tf_sound = cell(1,length(TF));
    tf_start_win = cell(1,length(TF));
    tf_sound_win = cell(1,length(TF));
    bl_l_trials = zeros(1,length(TF));
    for i = 1:length(TF) % trials in session
        tf_start{i} = TF(i).States.Start.sp(:,end - roi_wind+1:end); % window before sound
%         tf_start{i} = TF(i).States.Start.sp(:,1:roi_wind);
        stt = string(fieldnames(TF(i).States));
        fn = find(contains(stt , "Sound"));
        s_time = zeros(1,length(stt));
        for st = 1:length(stt)
            s_time(st) = floor(TF(i).States.(stt(st)).Tosca_Frames(1)*(twop_fr/tosca_fr));
        end
        s_time(fn) = s_time(fn) - 1;% possible rounding error fix
        lck = unique(ceil(find(TF(i).Lick)*(twop_fr/tosca_fr)));
        if ~isempty(lck)
            bl_licks = any(lck < s_time(fn));
            if bl_licks
                bl_l_trials(i) = 1;
                continue;
            end
        end
        if contains(TF(i).Group,'1')
            by_sound(1,i) = 1;
        else
            by_sound(2,i) = 1;
        end
        % if there are not enough frames in sound state, take from next
        % state
        if size(TF(i).States.(stt(fn)).sp,2) < roi_wind
            tf_nl = TF(i).States.(stt(fn+1)).sp;
            diff_frames = roi_wind - size(TF(i).States.(stt(fn)).sp,2);
            if size(TF(i).States.(stt(fn)).sp,2) + size(tf_nl,2) < roi_wind
                tf_sound{i} = nan;
                tf_start{i} = nan;
                continue;
            end
            tf_sound{i} = [TF(i).States.(stt(fn)).sp tf_nl(:,1:diff_frames)];
        else
            tf_sound{i} = TF(i).States.(stt(fn)).sp(:,1:roi_wind); % window after sound
        end
        % get frames of peak activity at each state
        [tf_start_win{i},tf_sound_win{i}] = peak_frames_ik(tf_start{i},tf_sound{i},sound_wind,method);
        %
%         figure(); plot(mean(TF(i).trial_sp));
%         xline(s_time(fn));
%         close;
        %
    end
    ind_bl_lick = logical(bl_l_trials);
    % get index of responsive cells
    % take only trials with quiet baseline (no licks)
    ind_trials = cellfun(@(x) ~isempty(x),tf_start_win);
    tf_start_win = tf_start_win(ind_trials);
    tf_sound_win = tf_sound_win(ind_trials);
    by_sound = by_sound(:,ind_trials);
    [responsive_cells,~,~] = respo_cells_ik(tf_start_win,tf_sound_win,split,by_sound);
    ret_responsive_cells = any(cell2mat(responsive_cells));

end