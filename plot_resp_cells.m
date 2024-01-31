function plot_resp_cells(resp_cells,sess_info)
    global num_grp
     % binned
    resp_bin = cell(1,length(num_grp));
    [m_bin,std_bin] = deal(zeros(length(num_grp),1));
    for b = 1:length(num_grp)
        for s_i = 1:length(resp_cells)
            if ismember(double(sess_info.Session(s_i)),num_grp(b,1):num_grp(b,2))
                resp_bin{b} = [resp_bin{b}; (sum(resp_cells{s_i})/length(resp_cells{s_i}))*100];
            end
        end
        m_bin(b) = mean(resp_bin{b});
        std_bin(b) = std(resp_bin{b})/sqrt(length(resp_bin{b}))
    end
    figure();
    bar(m_bin,'FaceColor','w')
    hold on
    errorbar(m_bin,std_bin,'LineStyle','none','Color','k');
    xlabel('Bin');
    ylabel('% responsive');
    box off
end