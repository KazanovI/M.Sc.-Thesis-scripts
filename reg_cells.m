function reg_cells(data,sess_info,dates_data,ind_cell_stats)
% data with original raw cell indices
% number of registered cells and their stat activity
    cd ('C:\onedriveitay\Master\Resnik lab\registrations\')
    mice = dir;
    mice_inds = find([mice.isdir]);
    m_names = string({mice(mice_inds(3:end)).name});
    all_mice_reg = arrayfun( @(x)struct('mouse',x,'cells',x,'dates',x,'target_column_ind',x),nan(1,length(m_names)));
    for m = 1:length(m_names)
        cd(['C:\onedriveitay\Master\Resnik lab\registrations\' char(m_names(m))])
        all_regis = dir('*.mat'); % load registerd files
        regis_name = string({all_regis.name});
        name_spl = regexp(regis_name,'_' ,'split');
        target_date = name_spl{1}(2); % extract dates and sort
        dates = cellfun(@(x) x(5),name_spl);
        all_dates = [target_date dates];
        dates_num = datenum(all_dates,'ddmmyy'); % to check chronological order
        [~,s_ind] = sort(dates_num);
        all_dates = all_dates(s_ind);
        target_ind = find(s_ind == 1);
        % get indices of cells
        curr_inds = cell(1,length(regis_name));
        for n = 1:length(regis_name)
            load(regis_name(n));
            curr_inds{n} = regi.rois.iscell_idcs;
        end
        max_ind_target = max(cellfun(@(x) max(x(:,1)),curr_inds)); % check maximal session num
        align_cells = nan(max_ind_target,length(curr_inds)+1);
        for i = 1:length(curr_inds)
            align_cells(curr_inds{i}(:,1),1) = curr_inds{i}(:,1);
            align_cells(curr_inds{i}(:,1),i+1) = curr_inds{i}(:,2);
        end
        align_cells =  align_cells(:,s_ind); % sort cells columns as well
        % if there are 3 or more registered cells *OR* 2 at first and last bin,
        % keep
        align_cell_pruned_ind = sum(~isnan(align_cells),2) > 2 |...
            (~isnan(align_cells(:,1)) & (~isnan(align_cells(:,end)) | ~isnan(align_cells(:,end-1))));
        align_cell_pruned = align_cells;
        align_cell_pruned(~align_cell_pruned_ind,:) = [];
        all_mice_reg(m).cells = align_cell_pruned; % cells 
        all_mice_reg(m).dates = all_dates;
        all_mice_reg(m).target_column_ind = target_ind;
        all_mice_reg(m).mouse = m_names(m);
    end

    sess_info = addvars(sess_info,dates_data','NewVariableNames','date');
    % add another field with reg cells data
    for i =  1:length(all_mice_reg)
        target_date = all_mice_reg(i).dates(all_mice_reg(i).target_column_ind);
        mice_data_sessions = cell(1,length(all_mice_reg(i).dates));
        name_sess = string();
        % per session, extract data
        for d = 1:length(all_mice_reg(i).dates) % sessions
            sess_date = all_mice_reg(i).dates(d); % current session date
            sess_ind = find(string(sess_info.date) == sess_date);
            if length(sess_ind) > 1
                % there are 2 recordings in the same date
                correct_mice = find(contains(sess_info.Mouse(sess_ind),all_mice_reg(i).mouse));
                sess_ind = sess_ind(correct_mice);
            end
            name_sess(d) = sess_info.Session(sess_ind);
            data_sess = data{sess_ind}; % currnet session data
            n_ind = find(all(isnan(data_sess),[2 3])); % nan cells
            n_ind_reg = ismember(all_mice_reg(i).cells(:,d),n_ind); % nan cells that were registerd
            all_mice_reg(i).cells(n_ind_reg,d) = nan; % update 
            sess_cells = all_mice_reg(i).cells(:,d); % current session cells
            dat_ind = find(~isnan(sess_cells));
            sess_cells(isnan(sess_cells)) = []; % remove nan cells if exist
            cells_data_sess = data_sess(sess_cells,:,:);
            cells_data_sess_full = [];
            cells_full = nan(length(all_mice_reg(i).cells(:,d)),size(data_sess,2),size(data_sess,3));
            cells_full(dat_ind,:,:) = cells_data_sess;
            cells_data_sess_full= cells_full;
            mice_data_sessions{d} = cells_data_sess_full; 
        end
        all_mice_reg(i).data_cells = mice_data_sessions;
        all_mice_reg(i).names = name_sess;
    end

    % load registered cells and check their change over bins
    figure('WindowState','maximized');
    cell_stat_type = cell(1,length(all_mice_reg)); % keep track on cells info(1-enh,2-sup,3-none)
    for m = 1:length(all_mice_reg) % per mouse
        mouse_name = all_mice_reg(m).mouse; % name
        mouse_inds_logical = strcmp(sess_info.Mouse,mouse_name);
        reg_sess_num = double(all_mice_reg(m).names); % session number
        reg_sess_ind = ismember(double(sess_info.Session(mouse_inds_logical)),reg_sess_num);
        mouse_ind = find(mouse_inds_logical);
        reg_cell_stats = ind_cell_stats(mouse_ind(reg_sess_ind));
        [~,srt] = sort(double(sess_info.Session(mouse_ind(reg_sess_ind))));
        reg_cell_stats = reg_cell_stats(srt);
        cell_stat_type{m} = nan(size(all_mice_reg(m).cells));
        for s = 1:length(reg_cell_stats) % per reg session
            non_nan_ind = ~isnan(all_mice_reg(m).cells(:,s));
            cells_ind = all_mice_reg(m).cells(non_nan_ind,s); % reg cells ind
            [row, col] = find(reg_cell_stats{s}(cells_ind,:)); % find type of modulation
            [~,s_idx] = sort(row);
            cell_stat_type{m}(non_nan_ind,s) = col(s_idx); % save cells info
        end
        cell_stat_type{m}(cell_stat_type{m} == 2) = -1; % binary coding
        cell_stat_type{m}(cell_stat_type{m} == 3) = 0;
        % per cell check first and last point to define cell's final status
        non_nan = cellfun(@(x) find(~isnan(x)), mat2cell(cell_stat_type{m},repelem(1,1,size(cell_stat_type{m},1))...
            ,size(cell_stat_type{m},2)), 'UniformOutput', false);
        min_max_ind{m} = cell2mat(cellfun(@(x) [min(x) max(x)], non_nan, 'UniformOutput', false));
        for c = 1:size(cell_stat_type{m},1)
            dat_cells(c,:) = cell_stat_type{m}(c,min_max_ind{m}(c,:));
        end
        % check type of change between first and last pnts
        % first, remove cells than were 0 at the end (no statistical modulation)
        end_zero_ind = dat_cells(:,2) == 0;
        dat_cells(end_zero_ind,:) = [];
        % second , cells that didnt change also remove
        same_ind = dat_cells(:,2) - dat_cells(:,1) == 0;
        same_cells = dat_cells(same_ind,:);
        same_enh = sum(same_cells(:,1) == 1);
        same_sup = size(same_cells,1) - same_enh;
        dat_cells(same_ind,:) = [];
        % last, check type of gain of function
        gain_enh = sum(dat_cells(:,2) == 1);
        gain_sup = sum(dat_cells(:,2) == -1);
        no_gain =   size(cell_stat_type{m},1) - (gain_enh + gain_sup);
        perc_gain = ([gain_enh same_enh gain_sup same_sup no_gain-size(same_cells,1)]/ size(cell_stat_type{m},1))*100;
        subplot(2,2,m) 
        lbls = ["Enhanced","No change(Enh)","Suppressed","No change(Sup)","None"];
        explode = [1 1 0 0 0];
        h = pie(perc_gain,explode,'%.2f%%');
        patchHand = findobj(h, 'Type', 'Patch'); 
        patchHand(1).FaceColor = 'b';
        patchHand(3).FaceColor = 'r';
        patchHand(5).FaceColor = [70,70,70]/255;
        patchHand(2).FaceColor = 'b'; patchHand(2).FaceAlpha = 0.2;
        patchHand(4).FaceColor = 'r'; patchHand(4).FaceAlpha = 0.2;
        % Create legend
        if m == 1
            legend(lbls, 'Location','northwestoutside');
        end
        title([all_mice_reg(m).mouse "    " ]);
    end
    sgtitle('Gain of function or no change, registered cells');

    % check all sessions of registered cells
    map = [181 60 60; 
           62 181 60;
           36 110 171]/255;
    figure();
    for m = 1:length(cell_stat_type)
        subplot(2,2,m)
        h = imagesc(cell_stat_type{m});
        no_data_clr = double(~isnan(cell_stat_type{m}));
        no_data_clr(no_data_clr==0) = 0.2;
        set(h, 'AlphaData', no_data_clr);
        colormap(map);
        colorbar('Ticks',[-1,0,1],...
         'TickLabels',{'Suppressed','None','Enhanced'});
     nms = all_mice_reg(m).names;
     xticks(1:length(nms));
     xticklabels(nms);
     xlabel("Session");
     ylabel("Cells");
     title(all_mice_reg(m).mouse);
    end
end