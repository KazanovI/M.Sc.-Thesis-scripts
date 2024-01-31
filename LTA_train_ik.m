function [data_both,inv_cells_sess,ret_trial_track,inv_trials_sess_both,...
    cnt_licks,ret_licks_time,ret_m_fr_sess,ind_resp] = LTA_train_ik(Tosca_F)
    % extract trials and check conditions
    [TF,trial_track] = extract_conds(Tosca_F,0);
    % Sound responsive cells
    [ind_resp] = sound_resp_cells(TF);
    % extrcat LTA 
    [sess_LTA_start,inv_trials_sess_start,cnt_licks_start] = LTA_start(TF);
    [sess_LTA_task,inv_trials_sess_task,cnt_licks_task,licks_time] = LTA_task(TF,1);   
    % indices of outlier cells
    [inv_cells_sess,m_fr_sess] = check_outlier_neurons([{sess_LTA_start}; {sess_LTA_task}]);
    % return data
    data_both = [{sess_LTA_start}; {sess_LTA_task}];
    inv_trials_sess_both = {[inv_trials_sess_start; inv_trials_sess_task]};
    cnt_licks = [{cnt_licks_start};{cnt_licks_task}];
    ret_licks_time = {licks_time};
    ret_trial_track = [{trial_track};{trial_track}]; % replicate for start/task tracking seperatly
    ret_m_fr_sess = {m_fr_sess};
end