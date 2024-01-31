function plot_start_task2_svm(svm_res_cond,varbs,varbs_check,cts,sv,shuf)
    cts_name = ["Enhanced","Suppressed","None"];
    num_iters = 5000;
    pv = permutationTest(svm_res_cond(:,1), svm_res_cond(:,2),num_iters);
    fprintf('start-task. p-value: %.03d, iters: %d',pv,num_iters);
    clrs = get_color(cts_name(cts));
    figure();
    distributionPlot(svm_res_cond(:,1),'widthDiv',[2 1],'histOri','left','color',clrs(1,:),'showMM',0,'xValues',1);
    hold on
    distributionPlot(svm_res_cond(:,2),'widthDiv',[2 2],'histOri','right','color',clrs(2,:),'showMM',0,'xValues',2);
    plot([1.2 1.8],mean(svm_res_cond),'xk','MarkerSize',12,'MarkerFaceColor','k');
    xticks(1:2)
    xticklabels(cts_name(cts));
%     title(['permutation test(1000 iter) p-value: ' num2str(pv)]);
    ylabel('prediction accuracy');
%     xlabel('Cell type');
    if ~shuf
        ylim([20 100]);
    else
        ylim([20 100]);
    end
    if sv
        save2 = ['C:\onedriveitay\Master\Resnik lab\results\SVM_figs'];
        if ~exist(save2, 'dir')
            mkdir(save2);
        end
        feat = char(strrep(strjoin(varbs(2,varbs_check))," ","_"));
        feat = string(feat(1:end-1));
        if shuf
            feat = append(feat,"shuff");
        end
        fig_name = append("new_train_test_",cts_name(cts(1)),"_",cts_name(cts(2)),"_",feat,".fig");
        saveas(gcf,fullfile(save2,fig_name));
        fig_name = append("new_train_test_",cts_name(cts(1)),"_",cts_name(cts(2)),"_",feat,".jpeg");
        saveas(gcf,fullfile(save2,fig_name));
    end
end