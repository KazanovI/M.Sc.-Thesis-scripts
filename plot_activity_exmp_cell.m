function plot_activity_exmp_cell(data,trial_track,info_sess,time,tst_cnd,sv)
    d = ['C:\onedriveitay\Master\Resnik lab\results\cell_activity\' time];
    mkdir(d);
    cd(d)
    if tst_cnd
        d = ['C:\onedriveitay\Master\Resnik lab\results\cell_activity\' time '\cnd_split'];
        mkdir(d);
        cd(d)
    end
    cnds_tst = [1 2; 3 4; 5 6];
    wnd = 28:35;
    nm_trials = 25;
    for s_i = 1:length(data)
        if info_sess.Phase(s_i) == "Test"
            if tst_cnd
%                 zcnd_sort = cell(1,3);
%                 for cnd = 1:3
%                     z_cnd = zscore(squeeze(mean(data{s_i}(:,:,ismember(trial_track{s_i},cnds_tst(cnd,:))),3)),[],2);
%                     [~,srt_i] = sort(mean(z_cnd(:,wnd),2),'descend');
%                     zcnd_sort{cnd} = z_cnd(srt_i,:);
%                 end
%                 zdata_sort = cell2mat(zcnd_sort');
%                 
% %                 title(['Mouse: ' char(info_sess.Mouse(s_i)), ', Session: ' char(info_sess.Session(s_i))]);
%                 yticks(cnds_ind);
%                 yticklabels(["Hit","Miss","FA"]);
%                 yline(yln(1:2),'LineWidth',5);
%                 box off;
%                 if sv
%                     saveas(gcf,[char(info_sess.Mouse(s_i)), '_Session_' char(info_sess.Session(s_i)),'.fig']);
%                     saveas(gcf,[char(info_sess.Mouse(s_i)), '_Session_' char(info_sess.Session(s_i)),'.jpg']);
%                 end
%                 close;
            else
            try
                zdata = mean(zscore(data{s_i}(:,:,1:nm_trials),[],[2 3]),3)';
                std_data = (std(zscore(data{s_i}(:,:,1:nm_trials),[],[2 3]),[],3)/sqrt(nm_trials))';
            catch
                zdata = mean(zscore(data{s_i},[],[2 3]),3)';
                std_data = (std(zscore(data{s_i},[],[2 3]),[],3)/sqrt(size(data{s_i},3)))';
            end
                [~,srt_i] = sort(mean(zdata(wnd,:),1),'descend');
                zdata_sort = zdata(:,srt_i);
                std_data_sort = std_data(:,srt_i);
                zdata_sort_sm = smoothdata(zdata_sort,1,'gaussian',3);
                std_data_sort_sm = smoothdata(std_data_sort,1,'gaussian',3);
                for c_i = 1:size(zdata_sort_sm,2)
                    % plot
                    figure('Visible','off');
%                     figure();
                    plot(zdata_sort_sm(:,c_i),'k');
                    hold on 
                    x2 = [1:length(zdata_sort_sm(:,c_i)), fliplr(1:length(zdata_sort_sm(:,c_i)))];
                    inBetween = [(zdata_sort_sm(:,c_i)-std_data_sort_sm(:,c_i)) ;...
                    flipud(zdata_sort_sm(:,c_i)+std_data_sort_sm(:,c_i))];
                    fill(x2, inBetween,'k','FaceAlpha',0.3,'EdgeAlpha',0.5);
%                     title([time,',Mouse: ' char(info_sess.Mouse(s_i)), ', Session: '...
%                     char(info_sess.Session(s_i)), ',Cell: ' num2str(c_i)]);
                    box off
                    ylim([-2.5 3.5]);
                    xticks(0:5:length(zdata_sort_sm(:,c_i)));
                    digits(2);
                    xticklabels(string(vpa(-30:5:15)/30));
                    xline(30,':k')
                    if sv
%                         saveas(gcf,[char(info_sess.Mouse(s_i)), '_Session_' char(info_sess.Session(s_i)), '_Cell_' num2str(c_i),'.fig']);
                        saveas(gcf,[char(info_sess.Mouse(s_i)), '_Session_' char(info_sess.Session(s_i)), '_Cell_' num2str(c_i),'.jpg']);
                    end
                    close;
                end
            end
        else
%             zdata = zscore(mean(data{s_i},3),[],2)';
%             std_data = (std(data{s_i},[],3)/sqrt(size(data{s_i},3)))';
            try
                zdata = mean(zscore(data{s_i}(:,:,1:nm_trials),[],[2 3]),3)';
                std_data = (std(zscore(data{s_i}(:,:,1:nm_trials),[],[2 3]),[],3)/sqrt(nm_trials))';
            catch
                zdata = mean(zscore(data{s_i},[],[2 3]),3)';
                std_data = (std(zscore(data{s_i},[],[2 3]),[],3)/sqrt(size(data{s_i},3)))';
            end
            [~,srt_i] = sort(mean(zdata(wnd,:),1),'descend');
            zdata_sort = zdata(:,srt_i);
            std_data_sort = std_data(:,srt_i);
            zdata_sort_sm = smoothdata(zdata_sort,1,'gaussian',3);
            std_data_sort_sm = smoothdata(std_data_sort,1,'gaussian',3);
            for c_i = 1:size(zdata_sort_sm,2)
                % plot
                figure('Visible','off');
%                 figure();
                plot(zdata_sort_sm(:,c_i),'k');
                hold on 
                x2 = [1:length(zdata_sort_sm(:,c_i)), fliplr(1:length(zdata_sort_sm(:,c_i)))];
                inBetween = [(zdata_sort_sm(:,c_i)-std_data_sort_sm(:,c_i)) ;...
                flipud(zdata_sort_sm(:,c_i)+std_data_sort_sm(:,c_i))];
                fill(x2, inBetween,'k','FaceAlpha',0.3,'EdgeAlpha',0.5);
%                 title([time,',Mouse: ' char(info_sess.Mouse(s_i)), ', Session: '...
%                 char(info_sess.Session(s_i)), ',Cell: ' num2str(c_i)]);
                box off
                ylim([-2.5 3.5]);
                xticks(0:5:length(zdata_sort_sm(:,c_i)));
                digits(2);
                xticklabels(string(vpa(-30:5:15)/30));
                xline(30,':k')
                if sv
%                     saveas(gcf,[char(info_sess.Mouse(s_i)), '_Session_' char(info_sess.Session(s_i)), '_Cell_' num2str(c_i),'.fig']);
                    saveas(gcf,[char(info_sess.Mouse(s_i)), '_Session_' char(info_sess.Session(s_i)), '_Cell_' num2str(c_i),'.jpg']);
                end
                close;
            end
        end
    end
end