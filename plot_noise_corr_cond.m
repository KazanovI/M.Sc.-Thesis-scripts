function [bin_cell_dat,bin_ind_cell] = plot_noise_corr_cond(data,sess_info,trial_track)   
    global num_grp
    cnds_tst = [1 2; 3 4; 5 6];
    window = [1:15;16:30]; % further/closer lick
%     window = [1:15;23:37]; % before/during lick
    bin_dat = cell(1,2);
    [bin_cell_dat,bin_ind_cell] = deal(cell(1,2));
    for ii = 1:size(window,1)
        bin_dat{ii} = cell(1,length(num_grp));
        [bin_cell_dat{ii},bin_ind_cell{ii}] = deal(cell(1,length(num_grp)));
        % all cells, by condition
        for b = 1:size(num_grp,1) % per bin
            bin_cor = []; % initialize
            cnt = 0;
            for s_i = 1:length(data)% session
                % check if in bin
                if ismember(double(sess_info.Session(s_i)),num_grp(b,1):num_grp(b,2))
                    dat_cat = data{s_i}(:,window(ii,:),:);
                    [bet_cor,bet_cor_all,ind_mat] = calc_bet_corr2(dat_cat);
                    % if test
                    if sess_info.Phase(s_i) == "Test"
                        tmp = [];
                        for cond = 1:3
                            ind = ismember(trial_track{s_i},cnds_tst(cond,:));
                            tmp(:,cond) = mean(bet_cor_all(:,ind),2,'omitnan');
                        end
                        bin_cor = [bin_cor ; tmp];
                    else
                        bin_cor = [bin_cor ; bet_cor];
                    end
                    cnt = cnt + 1;
                    bin_cell_dat{ii}{b}{cnt} = bet_cor;
                    bin_ind_cell{ii}{b}{cnt} = ind_mat;
                end
            end
            bin_dat{ii}{b} = bin_cor; 
        end
    end
    %ANOVA
    dat = cellfun(@(x) x(4:6),bin_dat,'UniformOutput',false);
    dat_all = [];
    time = [];
    cond = [];
    indx = [1 2 3;4 5 6];
    for i = 1:length(dat)
        dat_times{i} = cell2mat(dat{i}');
        time = [time; repelem(i,length(dat_times{i})*3)'];
        cond = [cond; repelem([1 2 3]',length(dat_times{i}))];
        tmp = dat_times{i};
        dat_all = [dat_all ; cat(1,tmp(:))];
    end
    figure();
    [p,tbl,stats] = anovan(dat_all,{time cond},'model',2,'varnames',{'time','cond'});
    [c,~,~,gnms] = multcompare(stats,'Dimension',[1 2],'CType','bonferroni');
    ns = string(gnms(c(c(:,6) > 0.05,1:2)))
    
    
    [mean_bin,std_bin] = deal(cell(1,2));
    for i = 1:length(bin_dat)
        tmp = cellfun(@(x) mean(x,'omitnan'),bin_dat{i},'UniformOutput',false);
        mean_bin{i} = cell2mat(tmp(4:6)');
        tmpp = cellfun(@(x) std(x,'omitnan')/sqrt(length(x)),bin_dat{i},'UniformOutput',false);
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
    ylim([0 0.03]);
    box off
    
    % cat bins
    [mean_bin,std_bin] = deal(cell(1,2));
    for i = 1:length(bin_dat)
        mean_bin{i} = mean(cell2mat(bin_dat{i}(4:6)'),'omitnan');
        std_bin{i} = std(cell2mat(bin_dat{i}(4:6)'),'omitnan')/sqrt(length(cell2mat(bin_dat{i}(4:6)')));
    end
    stds = cell2mat(std_bin');
    figure();
    dats = cell2mat(mean_bin');
    for i = 1:2
        br = bar(i,dats(i,:));
        if i == 1
            legend(br,["Hit","Miss","FA"],'AutoUpdate','off');
        end
        xends = [br.XEndPoints];
        yends = [br.YEndPoints];
        hold on
        errorbar(xends,yends,stds(i,:),'+k');
        for clr = 1:3
            br(clr).FaceColor = clrs(clr,:);
            if i == 2
                hatchfill2(br(clr),'cross','HatchAngle',45,'hatchcolor','k','HatchLineWidth',1);
            end
        end
    end
    xticks(1:2);
    xticklabels(categorical(["Further","Closer"]));
    box off
end