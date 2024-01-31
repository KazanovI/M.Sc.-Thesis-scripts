function num_sess_bin(indices,sess_info)
    global num_grp
    % binned
    num_bin = zeros(1,length(num_grp));
    for b = 1:length(num_grp)
        for s_i = 1:length(indices)
            if ismember(double(sess_info.Session(s_i)),num_grp(b,1):num_grp(b,2))
                num_bin(b) = num_bin(b) + 1;
            end
        end
    end
    figure();
    b = bar(num_bin,0.9,'k');
    xlabel("Bin");
    ylabel("Number of sessions");
    title("Number of sessions per bin");
    xtips1 = b.XEndPoints;
    ytips1 = b.YEndPoints;
    labels1 = string(b.YData);
    text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom')
    axis padded
    box off
end