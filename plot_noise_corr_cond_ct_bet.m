function plot_noise_corr_cond_ct_bet(data,sess_info,inds,trial_track)   
    global num_grp
    cnds_tst = [1 2; 3 4; 5 6];
    cts_names = ["Enh-Supp","Enh-None","Supp-None"];
    cts_ind = [1 2; 1 3; 2 3];
%     window = [1:15;16:30]; % before/during lick
    window = [1:15;23:37]; % before/during lick
    for cts = 1:size(inds{1},2)
        bin_dat = cell(1,2);
        for ii = 1:size(window,1) 
            bin_dat{ii} = cell(1,length(num_grp));
            % all cells, by condition
            for b = 1:size(num_grp,1) % per bin
                bin_cor = []; % initialize
                for s_i = 1:length(data)% session
                    % check if in bin
                    if ismember(double(sess_info.Session(s_i)),num_grp(b,1):num_grp(b,2))
                        ind_ct1 = logical(inds{s_i}(:,cts_ind(cts,1)));
                        dat_cat1 = data{s_i}(ind_ct1,window(ii,:),:);
                        ind_ct2 = logical(inds{s_i}(:,cts_ind(cts,2)));
                        dat_cat2 = data{s_i}(ind_ct2,window(ii,:),:);
                        if size(dat_cat1,1) <= 3 || size(dat_cat2,1) <= 3
                            continue;
                        else
                            [bet_cor,bet_cor_all] = calc_bet_corr2(dat_cat1,dat_cat2);
                        end
                        % if test
                        if sess_info.Phase(s_i) == "Test"
                            tmp = [];
                            for cond = 1:3
                                ind = ismember(trial_track{s_i},cnds_tst(cond,:));
                                both_cells = mean(bet_cor_all(:,ind,:),2,'omitnan');
                                tmp(:,cond) = both_cells(:);
                            end
                            bin_cor = [bin_cor ; tmp];
                        else
                            bin_cor = [bin_cor ; bet_cor];
                        end
                    end
                end
                bin_dat{ii}{b} = bin_cor; 
            end
        end
        [mean_bin,std_bin] = deal(cell(1,2));
        for i = 1:length(bin_dat)
            tmp = cellfun(@(x) mean(x),bin_dat{i},'UniformOutput',false);
            mean_bin{i} = cell2mat(tmp(4:6)');
            tmpp = cellfun(@(x) std(x)/sqrt(length(x)),bin_dat{i},'UniformOutput',false);
            std_bin{i} = cell2mat(tmpp(4:6)');
        end
        stds = cell2mat(std_bin);
        clrs = get_color("Hit","Miss","FA");
        clrs = repmat(clrs,2,1);
        figure();
        br = bar([1 2 3],cell2mat(mean_bin));
        for i = 1:length(br)
            xends = br(i).XEndPoints;
            yends = br(i).YEndPoints;
            hold on
            errorbar(xends,yends,stds(:,i),'+k');
            br(i).FaceColor = clrs(i,:);
            if i > 3
                hatchfill2(br(i),'cross','HatchAngle',45,'hatchcolor','k','HatchLineWidth',1);
            end
        end
        ylim([-0.03 0.025]);
        box off
        title(cts_names(cts));
    end
end