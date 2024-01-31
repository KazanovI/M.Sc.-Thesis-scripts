function [trials_ind] = clean_outlier_cov(dat_sess)
    % compute average of single-trial covariances
    % Without storing individual covmats!
    covave = zeros(size(dat_sess,1));
    for triali = 1:size(dat_sess, 3)
        covave = covave + cov( squeeze(dat_sess(:,:,triali))');
    end
    % divide by number of trials
    covave = covave / triali;
    % now loop through trials and compute the distance to the average
    covdist = zeros(size(dat_sess,3),1);
    for triali=1:size(dat_sess,3)
        thistrialcov = cov(squeeze(dat_sess(:,:,triali))');
        % compute Frobenius distance
        covdist(triali) = sqrt( sum(thistrialcov(:) .* covave(:)) );
    end
    % convert to z
    covdistz = (covdist-mean(covdist)) / std(covdist);
    %% pick a threshold and reject trials
    % threshold
    thresh = 2.3; % ~.01
    % identify trials that exceed the threshold
    trials_ind = abs(covdistz)>thresh;
end