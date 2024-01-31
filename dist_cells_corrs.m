function dist_cells_corrs(locs,data,inds,sess_info)
    global num_grp
    locs_bin = cell(1,length(num_grp));
    dist_bin = cell(1,length(num_grp));
    for b = 1:size(num_grp,1) % per bin
        cnt = 0;
        for s_i = 1:length(locs)% session
            % check if in bin
            if ismember(double(sess_info.Session(s_i)),num_grp(b,1):num_grp(b,2))
                cnt = cnt + 1;
                locs_bin{b}{cnt} = locs{s_i};
            end
        end
        for s_i = 1:length(data{2}{b})% session
            inds_sess = inds{2}{b}{s_i};
            dist = zeros(length(inds_sess),1);
            for cell_i = 1:length(dist) % pair-wise
                p_locs = locs_bin{b}{s_i}(inds_sess(cell_i,:),:);% pair locs
                dist(cell_i) = sqrt(diff(p_locs(:,1))^2 + diff(p_locs(:,2))^2);
            end
            dist_bin{b}{s_i} = dist;
        end
    end
    dist_all = cell2mat(cellfun(@(x) cell2mat(x'),dist_bin,'UniformOutput',false)');
    data_all = cell2mat(cellfun(@(x) cell2mat(x'),data{2},'UniformOutput',false)');
    [d_sorted,s_i] = sort(dist_all,'ascend');
    dat_sess_sort = data_all(s_i);
    [n,edgs] = histcounts(1:length(d_sorted),6);
    edgs(1) = 1; edgs(end) = length(d_sorted);
    dist_names = d_sorted(edgs);
    br = cell(1,length(n));
    for edg = 1:length(n)
        br{edg} = dat_sess_sort(edgs(edg):edgs(edg+1));
    end
    m_edg = cellfun(@(x) mean(x),br);
    std_edg = cellfun(@(x) std(x)/sqrt(length(x)),br);
    figure();
    bar(m_edg,'FaceColor',[0.6 0.6 0.6]);
    hold on 
    errorbar(m_edg,std_edg,'k');
    digits(2);
    xticklabels(string(vpa(dist_names)));
    xlabel('Distance between cells');
    ylabel('Noise corr');
    box off
    
    figure();
    scatter(d_sorted,abs(dat_sess_sort),1,'+k');
    box off
end