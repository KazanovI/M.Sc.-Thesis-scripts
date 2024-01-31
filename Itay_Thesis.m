%% Thesis Itay Kazanovich
% pipeline for all scripts and figures together
clear;clc;
%% sync and align data
sync = 0;
if sync
    Tosca_Imaging_Sync_ik;
end
%% extract Lick triggered average and general lick data
% load raw data
cd ('C:\onedriveitay\Master\Resnik lab\Sync\')
mice_synced = dir('**/*.mat');
% initialize
sessions = string({mice_synced.name}); % identifier
valid_sess_LTA = strings(length(sessions),4); % order
ind_inv_trials = cell(1,length(sessions)); % invalid trials, no lick 
inv_cells = cell(1,length(sessions)); % invalid cells, outlier FR
trial_track = cell(2,length(sessions)); % condition tracker
data_all = cell(2,length(sessions)); % LTA
licks_count = cell(2,length(sessions)); % binary, whether there were lick in trial 
lcks_rate = cell(1,length(sessions)); % frame of licks per trial
dates = NaT(1,length(sessions)); % date of session
FR = cell(1,length(sessions)); % mean firing rate per session
resp_cells = cell(1,length(sessions)); % indices of sound responsive cells
cells_locs = cell(1,length(sessions)); % locations of cells
cell_inh = cell(1,length(sessions)); % indices of inhibitory cells
% initialize globlas
global tosca_fr twop_fr wind_pre wind_post End_max_frames q_time mouse num_grp
wind_pre = 30; % frames
wind_post = 15;
q_time = 10; % num of frames w/o lick
End_max_frames = 350; % maximum frames of trial to take
tosca_fr = 500; % frame rate
twop_fr = 30; % fps
num_grp = [1 2 ; 3 4 ; 5 7 ; 8 10 ; 11 13 ; 14 17]; % bins
% load existing or pre-process new?
answer = inputdlg('load existing = 0, pre-process again = 1');
save_loc = 'C:\onedriveitay\Master\Resnik lab\processed results';
if answer{1} == '1'
     %save data after processing
    for ss = 1:length(sessions)
        keys = regexp(sessions{ss},'[_.]','split');
        mouse = char(keys(1));
        session = char(keys(3));
        %load session data
        cd (['C:\onedriveitay\Master\Resnik lab\Sync\' mouse '\Session ' num2str(session) '\']);
        load([mouse '_Session_' session]) % Tosca_F
        % return data by trial and invalid cells (outlier firing rate) and
        % trials
        if ~any(contains(string(Tosca_F.runs{1, 1}.trial{1, 1}.History),"NoLick")) % train sess
            [data_all(:,ss),inv_cells{ss},trial_track(:,ss),ind_inv_trials(ss),...
            licks_count(:,ss),lcks_rate(ss),FR(ss),resp_cells{ss}] = LTA_train_ik(Tosca_F);
            phase = 0; % train/test session
        else
            [data_all(:,ss),inv_cells{ss},trial_track(:,ss),ind_inv_trials(ss),...
            licks_count(:,ss),lcks_rate(ss),FR(ss),resp_cells{ss}] = LTA_test_ik(Tosca_F);
            phase = 1;
        end
        dates(ss) = datetime(Tosca_F.date,'InputFormat','dd/MM/yy','Format','ddMMyy');
        valid_sess_LTA(ss,1:3) = [string(mouse) ,string(session) , num2str(ss)];
        if phase
            valid_sess_LTA(ss,4) = "Test";
        else
            valid_sess_LTA(ss,4) = "Train";
        end
        cells_locs{ss} = Tosca_F.locations;
        cell_inh{ss} = Tosca_F.redCells;
    end
    valid_sess_LTA = array2table(valid_sess_LTA, 'VariableNames',{'Mouse','Session','Index','Phase'});
    dates = datetime(dates,'InputFormat','dd/MM/yy','Format','ddMMyy');
    %% pre-process: remove trials w/o lick, outlier cells, and outlier trials(by covariance matrix)
    % update trial track after cleaning
    data_all_clean = data_all; % initialize
    data_for_reg = data_all; % invalid cells set to nan
    trial_track_clean = trial_track;
    FR_clean = FR;
    for s_i = 1:length(data_all) % per session
        for phse = 1:size(data_all,1) % start/task
            % remove trials w/o lick
            ind_trials = ind_inv_trials{s_i}(phse,:);
            data_all_clean{phse,s_i}(:,:,logical(ind_trials)) = [];
            data_for_reg{phse,s_i}(:,:,logical(ind_trials)) = [];
            trial_track_clean{phse,s_i}(logical(ind_trials)) = [];
            if isempty(data_all_clean{phse,s_i})
                continue;
            end
            % remove outlier trial
            ind_trials_cov = clean_outlier_cov(data_all_clean{phse,s_i});
            f_ind = find(ind_trials_cov);
            for fg = 1:length(f_ind)
                figure();
                plot(data_all_clean{phse,s_i}(:,:,f_ind(fg))');
                close;
            end
            data_all_clean{phse,s_i}(:,:,ind_trials_cov) = [];
            data_for_reg{phse,s_i}(:,:,ind_trials_cov) = [];
            trial_track_clean{phse,s_i}(logical(ind_trials_cov)) = [];
            if isempty(data_all_clean{phse,s_i})
                continue;
            end
            % remove outlier cells
            ind_cells = inv_cells{s_i};
            data_all_clean{phse,s_i}(ind_cells,:,:) = [];
            data_for_reg{phse,s_i}(ind_cells,:,:) = nan;
        end
        FR_clean{s_i}(:,ind_cells) = [];
        resp_cells{s_i}(ind_cells) = [];
        cells_locs{s_i}(ind_cells,:) = [];
        cell_inh{s_i}(ind_cells,:) = [];
    end

    %% check stat cells
    % check significance of cells change in activity
    % cell can be enhanced, suppressed or none
    borders = [1 5; 30 34]; %equal lengths (frames)
    ind_cell_stats_start = cell_activity_stat(data_all_clean(1,:),borders);
    ind_cell_stats_task = cell_activity_stat(data_all_clean(2,:),borders);
    % reg cells
    ind_cell_stats_reg = cell_activity_stat(data_for_reg(2,:),borders);
    % cells that are significant in start & task
    ind_cell_stats_shared = shared_cells(ind_cell_stats_start,ind_cell_stats_task);
    % cells that are sound and lick responsive
    ind_sound_lick_shared = shared_sound_lick(resp_cells,ind_cell_stats_task);
    % in which group are inhibitory cells
    ind_inh_stat_shared = shared_inh_stat(cell_inh,ind_cell_stats_task);
    % save
    save(append(save_loc,'\full_data'),'data_all_clean','valid_sess_LTA','trial_track','licks_count','lcks_rate',...
    'ind_cell_stats_task','ind_cell_stats_start','ind_cell_stats_reg','FR_clean','ind_cell_stats_reg',...
    'ind_cell_stats_shared','ind_sound_lick_shared','ind_inh_stat_shared','resp_cells','cells_locs');
elseif answer{1} == '0'
    % load data after processing
    load(append(save_loc,'\full_data'));
end
%% figures
%% figure 1 - 2p explanation,task explanation,descriptives(FR, num of cells per condition,
% licks), cells example
% num sessions per bin
num_sess_bin(ind_cell_stats_task,valid_sess_LTA);
% num cells total
num_cells_total = sum(cellfun(@(x) length(x),ind_cell_stats_task));
% pe cell type
num_cells_ct = sum(cell2mat(cellfun(@(x) sum(x),ind_cell_stats_task,'UniformOutput',false)'));
% lick data, *many plots* (per mice/session)
[lick_prob] = plot_lick_data(licks_count,lcks_rate,valid_sess_LTA,trial_track(2,:));
% lick rate bar over sessions
plot_lck_rate(lcks_rate,trial_track(2,:));
% number of trials in each condition
plot_trial_num(data_all_clean(1,:),trial_track_clean(1,:),valid_sess_LTA); % start
plot_trial_num(data_all_clean(2,:),trial_track_clean(2,:),valid_sess_LTA); % task
% plot mean fr of all session(lick)
plot_fr(FR_clean,valid_sess_LTA);
% plot(and save) population activity examples
test_cnd = 0; % split plot to conditions?
save = 0; % save figs to folder(path in function)?
plot_activity_exmp(data_all_clean(1,:),trial_track_clean(1,:),valid_sess_LTA,'start',test_cnd,save); % start
plot_activity_exmp(data_all_clean(2,:),trial_track_clean(2,:),valid_sess_LTA,'task',test_cnd,save); % task
plot_activity_exmp(data_all_clean(1,:),trial_track_clean(2,:),valid_sess_LTA,'task',test_cnd,save); % mix
% plot single cell activity examples
test_cnd = 0; % split plot to conditions?
save = 0; % save figs to folder(path in function)?
plot_activity_exmp_cell(data_all_clean(1,:),trial_track_clean(1,:),valid_sess_LTA,'start',test_cnd,save); % start
plot_activity_exmp_cell(data_all_clean(2,:),trial_track_clean(2,:),valid_sess_LTA,'task',test_cnd,save); % start
plot_resp_cells(resp_cells,valid_sess_LTA);
plot_lick_snd_shr(ind_sound_lick_shared);
%% figure 2 - stat cells(task+start) FR and count per bin,registerd cells
% registered cells data
reg_cells(data_for_reg(2,:),valid_sess_LTA,dates,ind_cell_stats_reg);
% plot mean fr of all session(start->end,not lick),cell types in single plot 
% plot_fr_ct(FR_clean,valid_sess_LTA,ind_cell_stats_start); % start
% plot_fr_ct(FR_clean,valid_sess_LTA,ind_cell_stats_task); % task
% plot mean fr of lick window
[fr_start_lick_bin] = plot_fr_lick(data_all_clean(1,:),valid_sess_LTA,ind_cell_stats_start,"start");% start
[fr_task_lick_bin] = plot_fr_lick(data_all_clean(2,:),valid_sess_LTA,ind_cell_stats_task,"task");% task
[fr_mix_lick_bin] = plot_fr_lick(data_all_clean(1,:),valid_sess_LTA,ind_cell_stats_task,"mix");% mix
% number & percentage of cells per cell type
[perc_cells_start] = plot_stat_cells(ind_cell_stats_start,valid_sess_LTA,"Start"); % start
[perc_cells_task] = plot_stat_cells(ind_cell_stats_task,valid_sess_LTA,"Task"); % task
plot_stat_cells(ind_cell_stats_shared,valid_sess_LTA,"Task"); % shared
% activity of stat cells
plot_stat_activity(data_all_clean(1,:),valid_sess_LTA,ind_cell_stats_start,"Start"); % start
act_task = plot_stat_activity(data_all_clean(2,:),valid_sess_LTA,ind_cell_stats_task,"Task"); % task
act_mix = plot_stat_activity(data_all_clean(1,:),valid_sess_LTA,ind_cell_stats_task,"mix"); % mix
plot_enh_peak(act_task,act_mix);
% plot single cell activity examples
test_cnd = 0; % split plot to conditions?
save = 1; % save figs to folder(path in function)?
plot_activity_exmp_cell_stat(data_all_clean(1,:),trial_track_clean(1,:),ind_cell_stats_task,...
    valid_sess_LTA,'mix',test_cnd,save); % mix (cells task, activity start)
plot_activity_exmp_cell_stat(data_all_clean(2,:),trial_track_clean(2,:),ind_cell_stats_task,...
    valid_sess_LTA,'task',test_cnd,save); % task

% stat cells distance 
stat_cell_dist(ind_cell_stats_task,cells_locs);
%% figure 3 - noise correlations within group(task+start+mix)
% noise correlations between group(task+start+mix)

% by condition, 2 windows before lick
[dat_bins,ind_cell_bins] = plot_noise_corr_cond(data_all_clean(2,:),valid_sess_LTA,trial_track_clean(2,:)); % task
plot_noise_corr_cond(data_all_clean(1,:),valid_sess_LTA,trial_track_clean(1,:)); % start
% by condition, before lick in moving window
% wind_size = 10;
% plot_noise_corr_cond_steps(data_all_clean(2,:),valid_sess_LTA,trial_track_clean(2,:),wind_size);
% by condition and cell type, before lick
plot_noise_corr_cond_ct(data_all_clean(2,:),valid_sess_LTA,ind_cell_stats_task,trial_track_clean(2,:));% task
% by condition and between cell types, before lick
% plot_noise_corr_cond_ct_bet(data_all_clean(2,:),valid_sess_LTA,ind_cell_stats_task,trial_track_clean(2,:));
% distance between cells and corrs
dist_cells_corrs(cells_locs,dat_bins,ind_cell_bins,valid_sess_LTA);


% within
plot_noise_corr(data_all_clean(1,:),valid_sess_LTA,ind_cell_stats_start,trial_track_clean(1,:),"start");
plot_noise_corr(data_all_clean(2,:),valid_sess_LTA,ind_cell_stats_task,trial_track_clean(2,:),"task");
plot_noise_corr(data_all_clean(1,:),valid_sess_LTA,ind_cell_stats_task,trial_track_clean(1,:),"mix");
% between
% plot_noise_corr_bet(data_all_clean(1,:),valid_sess_LTA,ind_cell_stats_start,"start");
% plot_noise_corr_bet(data_all_clean(2,:),valid_sess_LTA,ind_cell_stats_task,"task");
% plot_noise_corr_bet(data_all_clean(1,:),valid_sess_LTA,ind_cell_stats_task,"mix");

%% figure 4 - PCA(variance explained,num pcs) + trajectories(task+start) + distances
% by cell type
save_figs = 0;
cond_comp = [1 2]; % conditions to compare (hit miss fa)
plot_traj(data_all_clean(1,:),valid_sess_LTA,ind_cell_stats_start,trial_track_clean(1,:),"start",cond_comp,save_figs);
[dists_task] = plot_traj(data_all_clean(2,:),valid_sess_LTA,ind_cell_stats_task,trial_track_clean(2,:),"task",cond_comp,save_figs);
[dists_mix] = plot_traj(data_all_clean(1,:),valid_sess_LTA,ind_cell_stats_task,trial_track_clean(1,:),"mix",cond_comp,save_figs);
plot_diff_dists(dists_task,dists_mix,cond_comp);
% all cells
save_figs = 0;
cond_comp = [1 3]; % conditions to compare (hit miss fa)
var_exp_start = plot_traj_all(data_all_clean(1,:),valid_sess_LTA,trial_track_clean(1,:),"start",cond_comp,save_figs);
var_exp_task = plot_traj_all(data_all_clean(2,:),valid_sess_LTA,trial_track_clean(2,:),"task",cond_comp,save_figs);
plot_perc_var_exp(var_exp_start,var_exp_task);
% random sampling with return from trials, permutation test
var_exp_task2 = plot_traj_all2(data_all_clean(2,:),valid_sess_LTA,trial_track_clean(2,:),"task",cond_comp,save_figs);

%% figure 5 - SVM
% SVM1(predict if lick in start or task, by activity of cells after pca)
% per session, mean over frames of lick
[pa,pa_shuf] = svm_start_task(data_all_clean);
sv = 0;
plot_start_task_svm(pa,pa_shuf,valid_sess_LTA,sv);
% SVM1(predict if lick in start or task, by activity, licks, and perc cells)
% per session. activity zscored per cell type and averaged to single value
% per session
for i = 1:length(perc_cells_task)
    dat1 = mat2cell(perc_cells_task{i}',[1 1 1])';
    perc_task(:,i) = dat1;
    dat2 = mat2cell(perc_cells_start{i}',[1 1 1])';
    perc_start(:,i) = dat2;
end

% varbs = {1 2 3; 'licks ','perc cells ','mean FR '};
% varbs_check = [1 2 3];
% cts = [1 2];
% shuf = 0;
% sv = 0;
% [ppa2] = svm_start_task2(lick_prob,perc_start,perc_task,fr_start_lick_bin,fr_task_lick_bin,varbs,varbs_check,cts,shuf);
% plot_start_task2_svm(ppa2,varbs,varbs_check,cts,sv,shuf);

varbs = {1 2 3; 'licks ','perc cells ','Activity '};
varbs_check = [ 3];
cts = [1 3];
shuf = 0;
sv = 0;
[ppa2] = svm_start_task3(lick_prob,perc_start,perc_task,data_all_clean,[ind_cell_stats_start;ind_cell_stats_task],varbs,varbs_check,cts,shuf);
plot_start_task2_svm(ppa2,varbs,varbs_check,cts,sv,shuf);


% SVM2(predict if lick in train or test)
varbs = {1 2 3; 'licks ','perc cells ','Activity '};
varbs_check = [1 2 3];
cts = [1 2];
sv = 1;
shuf = 1;
[pa2] = svm_train_test(lick_prob,perc_task,data_all_clean(2,:),[ind_cell_stats_start;ind_cell_stats_task],varbs,varbs_check,cts,shuf);
plot_train_test_svm(pa2,varbs,varbs_check,cts,sv);


