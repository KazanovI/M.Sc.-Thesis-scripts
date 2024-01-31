function plot_fr(FR,sess_info)
    global num_grp
    % all cells over sessions
    clrs = [0.6 0.6 0.6; 0 0 0];
    mean_hist = zeros(size(FR{1},1),1);
    hist_ti = cell(1,size(FR{1},1));
    figure();
    for ti = 1:size(FR{1},1)
        dat_ti = cellfun(@(x) x(ti,:),FR,'UniformOutput',false);
        hist_ti{ti} = cell2mat(dat_ti);
        mean_hist(ti) = mean(hist_ti{ti});
        histogram(hist_ti{ti},100,'FaceColor',clrs(ti,:),'FaceAlpha',0.3,...
            'Normalization','probability'); %,'Orientation','horizontal'
%         title('All cells');
%         xlabel('FR(Hz)');
%         ylabel('Number of cells');
        hold on
    end
    legend("Start","Task");
    box off
    [h1,p1] = kstest2(hist_ti{1},hist_ti{2});
    % binned
    fr_bin = cell(1,length(num_grp));
    [m_bin,std_bin] = deal(zeros(2,length(num_grp)));
    for b = 1:length(num_grp)
        for s_i = 1:length(FR)
            if ismember(double(sess_info.Session(s_i)),num_grp(b,1):num_grp(b,2))
                fr_bin{b} = [fr_bin{b} FR{s_i}];
            end
        end
        m_bin(:,b) = mean(fr_bin{b},2);
        std_bin(:,b) = std(fr_bin{b},[],2)./sqrt(length(fr_bin{b}));
    end
    
    % plot per bin
    figure();
    for ti = 1:size(m_bin,1)
        fi(ti) = plot(m_bin(ti,:),'Color',clrs(ti,:),'LineWidth',2);
        hold on 
        x2 = [1:length(m_bin(ti,:)), fliplr(1:length(m_bin(ti,:)))];
        inBetween = [(m_bin(ti,:)-std_bin(ti,:)) fliplr(m_bin(ti,:)+std_bin(ti,:))]; % vector,same length as x2
        fill(x2, inBetween,clrs(ti,:),'FaceAlpha',0.3,'EdgeAlpha',0.5);
    end
%     title("Mean FR per bin");
%     xlabel("Bin");
%     ylabel("Mean FR(Hz)");
%     xlim([0 length(num_grp)+1]);
    xticks(1:length(num_grp));
    legend(fi,"Pre-task","Task");
    axis padded
    box off
end