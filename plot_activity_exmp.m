function plot_activity_exmp(data,trial_track,info_sess,time,tst_cnd,sv)
    d = ['C:\onedriveitay\Master\Resnik lab\results\activitymap\' time];
    mkdir(d);
    cd(d)
    if tst_cnd
        d = ['C:\onedriveitay\Master\Resnik lab\results\activitymap\' time '\cnd_split'];
        mkdir(d);
        cd(d)
    end
    cnds_tst = [1 2; 3 4; 5 6];
    wnd = 28:35;
    map = [0.0392 0.4314 0.7412
        0.3529 0.5882 0.8902
        0.6275 0.7490 0.8784
        1 1 1
        1.0000 0.8510 0.7529
        0.9098 0.2706 0.2706
        0.5647 0.2157 0.2863];
    for s_i = 1:length(data)
        if info_sess.Phase(s_i) == "Test"
            if tst_cnd
                zcnd_sort = cell(1,3);
                [cnds_loc,cnds_ind,yln] = deal(zeros(1,3));
                for cnd = 1:3
                    z_cnd = zscore(squeeze(mean(data{s_i}(:,:,ismember(trial_track{s_i},cnds_tst(cnd,:))),3)),[],2);
                    [~,srt_i] = sort(mean(z_cnd(:,wnd),2),'descend');
                    zcnd_sort{cnd} = z_cnd(srt_i,:);
                    if cnd == 1
                        cnds_loc(cnd) = size(zcnd_sort{cnd},1);
                        cnds_ind(cnd) = round(median([0 cnds_loc(cnd)]));
                        yln(cnd) = size(zcnd_sort{cnd},1);
                    else
                        cnds_loc(cnd) = size(zcnd_sort{cnd},1);
                        cnds_ind(cnd) = round(median([sum(cnds_loc(1:cnd-1)) sum(cnds_loc(1:cnd-1))+cnds_loc(cnd)]));
                        yln(cnd) = sum([yln(cnd-1) cnds_loc(cnd)]);
                    end
                end
                zdata_sort = cell2mat(zcnd_sort');
                imagesc(zdata_sort);
                colormap(map)
                caxis([-2.5 2.5]);
                colorbar;
                title(['Mouse: ' char(info_sess.Mouse(s_i)), ', Session: ' char(info_sess.Session(s_i))]);
                yticks(cnds_ind);
                yticklabels(["Hit","Miss","FA"]);
                digits(2);
                xticks(0:5:40);
                xticklabels(string(vpa([-30:5:15]/30)));
                yline(yln(1:2),'LineWidth',5);
                box off;
                if sv
                    saveas(gcf,[char(info_sess.Mouse(s_i)), '_Session_' char(info_sess.Session(s_i)),'.fig']);
                    saveas(gcf,[char(info_sess.Mouse(s_i)), '_Session_' char(info_sess.Session(s_i)),'.jpg']);
                end
                close;
            else
                zdata = zscore(mean(data{s_i},3),[],2)';
                [~,srt_i] = sort(mean(zdata(wnd,:),1),'descend');
                zdata_sort = zdata(:,srt_i);
                imagesc(zdata_sort');
                colormap(map)
                caxis([-2.5 2.5]);
                colorbar;
%                 title(['Mouse: ' char(info_sess.Mouse(s_i)), ', Session: ' char(info_sess.Session(s_i))]);
                digits(2);
                xticks(0:5:40);
                xticklabels(string(vpa([-30:5:15]/30)));
                box off
                if sv
                    saveas(gcf,[char(info_sess.Mouse(s_i)), '_Session_' char(info_sess.Session(s_i)),'.fig']);
                    saveas(gcf,[char(info_sess.Mouse(s_i)), '_Session_' char(info_sess.Session(s_i)),'.jpg']);
                end
                close;
            end
        else
            zdata = zscore(mean(data{s_i},3),[],2)';
            [~,srt_i] = sort(mean(zdata(wnd,:),1),'descend');
            zdata_sort = zdata(:,srt_i);
            imagesc(zdata_sort');
            colormap(map)
            caxis([-2.5 2.5]);
            colorbar;
%             title(['Mouse: ' char(info_sess.Mouse(s_i)), ', Session: ' char(info_sess.Session(s_i))]);
            digits(2);
            xticks(0:5:40);
            xticklabels(string(vpa([-30:5:15]/30)));
            box off
            if sv
                saveas(gcf,[char(info_sess.Mouse(s_i)), '_Session_' char(info_sess.Session(s_i)),'.fig']);
                saveas(gcf,[char(info_sess.Mouse(s_i)), '_Session_' char(info_sess.Session(s_i)),'.jpg']);
            end
            close;
        end
    end
end