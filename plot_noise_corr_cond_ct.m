function plot_noise_corr_cond_ct(data,sess_info,inds,trial_track)   
    global num_grp
    cnds_tst = [1 2; 3 4; 5 6];
%     cts_names = ["Enhanced","Suppressed","None"];
    window = [1:15;16:30]; % further/closer lick
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
                        ind_ct = logical(inds{s_i}(:,cts));
%                         %% look at data
%                         temp = zscore(mean(data{s_i}(ind_ct,:,:),3),[],2)';
%                         [ri,ci] = find(tril(abs(corrcoef(temp)),-1) < 0.1 &...
%                             tril(abs(corrcoef(temp)),-1) ~= 0);
%                         [ri,ci] = find(tril(abs(corrcoef(temp)),-1) > 0.6);
%                         if ~isempty(ri)
%                             for cori = 1:length(ri)
%                                 figure();
%                                 plot(temp(:,ri(cori)),'Color',[0.2667 0.4667 0.8078],'LineWidth',1.5);
%                                 hold on
%                                 plot(temp(:,ci(cori)),'Color',[0.9451 0.3529 0.3490],'LineWidth',1.5);
%                                 box off
%                                 xticks(0:5:45);
%                                 digits(2);
%                                 xticklabels(string(vpa(-30:5:15)/30));
%                                 ylim([-3 4]);
%                                 close;
%                             end
%                         end
%                         %%
                        dat_cat = data{s_i}(ind_ct,window(ii,:),:);
                        if size(dat_cat,1) <= 3
                            continue;
                        else
                            [bet_cor,bet_cor_all] = calc_bet_corr2(dat_cat);
                        end
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
                    end
                end
                bin_dat{ii}{b} = bin_cor; 
            end
        end
        % ANOVA2
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
        [p,tbl,stats] = anovan(dat_all,{time cond},'model',2,'varnames',{'time','cond'});
        [c,~,~,gnms] = multcompare(stats,'Dimension',[1 2],'CType','bonferroni');
        ns = string(gnms(c(c(:,6) > 0.05,1:2)));
%         for n = 1:size(ns,1)
%             fprintf('NS %s: %s and %s \n',string(cts),ns(n,1),ns(n,2));
%         end
        % plot
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
        ylim([0 0.1]);
        box off
        
        % cat bins
        [mean_bin,std_bin] = deal(cell(1,2));
        for i = 1:length(bin_dat)
            mean_bin{i} = mean(cell2mat(bin_dat{i}(4:6)'));
            std_bin{i} = std(cell2mat(bin_dat{i}(4:6)'))/sqrt(length(cell2mat(bin_dat{i}(4:6)')));
        end
        stds = cell2mat(std_bin');
        figure();
        dats = cell2mat(mean_bin');
        for i = 1:2
            br = bar(i,dats(i,:));
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
        ylim([0 0.07]);
        box off
    end
end