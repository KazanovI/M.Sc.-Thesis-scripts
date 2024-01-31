function [inv_cells_both,ret_m_fr] = check_outlier_neurons(data)
% check_outlier_neurons calculate mean firing rate and search outliers
% outlier firing rate defined in function
    global twop_fr
    % remove low spiking neurons
    sp_th = [1 100]; % spikes per second , to remove outlier spiking neurons
    inv_cells = logical(zeros(length(data),size(data{1},1)));
    for t = 1:length(data)
        spk = reshape(data{t},size(data{t},1),[]);
        sess_len = length(spk)/twop_fr; % seconds
        sess_len_f = floor(sess_len); 
        ses_rem = length(spk) - (sess_len_f*twop_fr); %remove some frames 
        % so its possible to reshape (less than a second removed from the end)
        sess_len_round = spk(:,1:end - ses_rem );
        fr = squeeze(sum(reshape(sess_len_round,size(sess_len_round,1),twop_fr,[]),2,'omitnan')); % mean firing rate
        m_fr = mean(fr(:,any(fr)),2,'omitnan'); % mean firing rate
        std_fr = std(fr(:,any(fr)),[],2)/size(fr(:,any(fr)),2);
        outlier_ind = (m_fr < sp_th(1) | m_fr > sp_th(2)) | (std_fr > 3); %ind of low spiking neurons and high std over trials
        outlier_ind = isnan(m_fr) | outlier_ind;
        inv_cells(t,:) = outlier_ind;
        ret_m_fr(t,:) = m_fr;
    end
    inv_cells_both = any(inv_cells);
end