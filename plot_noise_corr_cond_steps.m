function plot_noise_corr_cond_steps(data,sess_info,trial_track,wind_size)   
    global num_grp
    cnds_tst = [1 2; 3 4; 5 6];
    window = [1:wind_size]; % before/during lick
    step_size = 1;
    num_steps = (30-wind_size)/step_size;
    bin_steps = cell(1,num_steps);
    for ii = 1:num_steps
        bin_steps{ii} = cell(1,length(num_grp)/2);
        % all cells, by condition
        for b = 4:size(num_grp,1) % per bin
            bin_cor = []; % initialize
            for s_i = 1:length(data)% session
                % check if in bin
                if ismember(double(sess_info.Session(s_i)),num_grp(b,1):num_grp(b,2))
                    dat_cat = data{s_i}(:,window,:);
                    [bet_cor,bet_cor_all] = calc_bet_corr2(dat_cat);
                    % if test
                    if sess_info.Phase(s_i) == "Test"
                        tmp = [];
                        for cond = 1:3
                            ind = ismember(trial_track{s_i},cnds_tst(cond,:));
                            tmp(:,cond) = mean(bet_cor_all(:,ind),2,'omitnan');
                        end
                        bin_cor = [bin_cor ; tmp];
                    end
                end
            end
            bin_steps{ii}{b-3} = bin_cor; 
        end
        window = window + step_size; % step
    end
    
    [mean_bin,std_bin] = deal(cell(1,length(bin_steps)));
    for i = 1:length(bin_steps)
        tmp = cellfun(@(x) mean(x,'omitnan'),bin_steps{i},'UniformOutput',false);
        mean_bin{i} = cell2mat(tmp');
        tmpp = cellfun(@(x) std(x,'omitnan')/sqrt(length(x)),bin_steps{i},'UniformOutput',false);
        std_bin{i} = cell2mat(tmpp');
    end
    [cond_mean,cond_std] = deal(cell(1,length(bin_steps{1})));
    for i = 1:length(mean_bin{1}) % per cond and bin
        cond_mean{i} = cell2mat(cellfun(@(x) x(:,i),mean_bin,'UniformOutput',false));
        cond_std{i} = cell2mat(cellfun(@(x) x(:,i),std_bin,'UniformOutput',false));
    end
    tit_name = ["Hit","Miss","FA"];
    clrs = get_color(tit_name);
    figure();
    for i = 1:length(cond_mean) % condition
        ax(i) = subplot(length(cond_mean),1,i);
        pp = plot(cond_mean{i}','color',clrs(i,:),'LineWidth',2);
        pp(1).LineStyle = '-';
        pp(2).LineStyle = ':';
        pp(3).LineStyle = '-.';
        legend(pp,["4","5","6"]);
        title(tit_name(i));
        box off
    end
    ylim([0 0.02]);
    linkaxes(ax,'y');
    xlabel("Steps");
    ylabel("noise corr");
end