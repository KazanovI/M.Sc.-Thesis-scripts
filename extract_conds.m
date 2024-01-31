function [TF,trial_track] = extract_conds(Tosca_F,stage)
% extract_conds check indices of each condition
% Tosca_F - session data
% stage - 0 = Train, 1 = Test
    TF = [];
    for i = 1:length(Tosca_F.runs)
        is_s = cellfun(@(x) isstruct(x), Tosca_F.runs{1,i}.trial ); % check if struct
        if any(is_s==0)
            remve = (is_s==0);
            TF = [TF Tosca_F.runs{1,i}.trial{1,~remve}]; % remaining structs
        else
            TF = [TF Tosca_F.runs{1,i}.trial{1,:}]; % remaining structs
        end
    end
    if stage
        % conditions
        TF_h = {TF(:).History}; % extract history
        hit_trials = cellfun(@(x) sum(contains(x,'Lick'))==3, TF_h, 'UniformOutput', false);
        miss_trials  = cellfun(@(x) sum(contains(x,'Lick'))==2 & sum(contains(x,'TO'))==0 , TF_h, 'UniformOutput', false);
        hit_trials = [hit_trials{:}];
        miss_trials = [miss_trials{:}];
        %by sound
        TF_g = {TF(:).Group}; % extract sound
        one = cellfun(@(x) x == "sound1", TF_g, 'UniformOutput', false);
        one = [one{:}];
        two = ~one;
        % include FA
        fa_trials = cellfun(@(x) sum(contains(x,'TO'))==2, TF_h, 'UniformOutput', false);
        fa_trials = [fa_trials{:}];
        gf = [hit_trials;miss_trials;fa_trials;one;two];
        gf_group = [gf(1,:) == 1 & gf(4,:) == 1; gf(1,:) == 1 & gf(5,:) == 1; ...
            gf(2,:) == 1 & gf(4,:) == 1;gf(2,:) == 1 & gf(5,:) == 1; ...
            gf(3,:) == 1 & gf(4,:) == 1; gf(3,:) == 1 & gf(5,:) == 1]; % hit1;hit2;miss1;miss2;fa1;fa2
        [trial_track,~] = find(gf_group);
        if ~isrow(trial_track)
            trial_track = trial_track';
        end
    else
        for t = 1:length(TF)
            % get group
            g = TF(t).Group;
            cond_type = regexp(g,'\d','match');
            g_num = string(cond_type); % 1 | 2
            trial_track(t) = str2double(cell2mat(cond_type));
        end
    end
end