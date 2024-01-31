function [svm_pa,svm_pa_shuf] = svm_start_task(data_both)
    %% cat *all* cell types, reduce dimentionality of zscored data after mean, and then SVM
    m_wind = 30:35; 
    [svm_res_time,svm_res_time_shuf] = deal(cell(1,length(data_both)));
    [pred_acc,ftr_num,obs_num] = deal([]);

    for si = 1:size(data_both,2)
        nn_inds = zeros(2,3);
        dat = data_both(:,si);
%         zscore by mean over trials baseline activity
%         for time = 1:length(dat)
%             z_dat = dat{time};
%             mn = mean(z_dat(:,1:5,:),[2 3]);
%             sd = std(z_dat(:,1:5,:),[],[2 3])/sqrt(size(z_dat,1));
%             ztemp = zeros(size(z_dat));
%             for tr = 1:size(z_dat,3)
%                 ztemp(:,:,tr) = (z_dat(:,:,tr) - mn)./sd;
%             end
%             dat{time} = ztemp;
%         end
        % mean
        dat_m = cellfun(@(x) squeeze(mean(x(:,m_wind,:),2)),dat,'UniformOutput',false) ;
        % make a table with labels
        dataset = array2table([dat_m{1}';dat_m{2}']);
        label0 = repmat("Task",size(dat_m{2},2),1);
        label1 = repmat("Start",size(dat_m{1},2),1);
        dataset = addvars(dataset, [label0;label1],...
            'NewVariableNames','label');
        labels = dataset(:,end);
        t = tabulate(dataset.label);
        uniqueLabels = string(t(:,1));
        labelCounts = cell2mat(t(:,2));
        % to which group synth data is needed
        [~,mi] = min(labelCounts);
        % how much to add
        dif_labels = abs(diff(labelCounts))-1;
        if mi==1
            num2Add = [dif_labels,0];
        else
            num2Add = [0,dif_labels];
        end
        % algorithm params
        algorithm = "Safe-level SMOTE";
        k = 3;
%         rng(1);
        newdata = table;
        % for each class
        for ii=1:length(uniqueLabels)
            switch algorithm
                case "SMOTE"
                    [tmp,visdata] = mySMOTE(dataset,uniqueLabels(ii),num2Add(ii),...
                        "NumNeighbors",k, "Standardize", true);
                case "ADASYN"
                    [tmp,visdata]  = myADASYN(dataset,uniqueLabels(ii),num2Add(ii),...
                        "NumNeighbors",k, "Standardize", true);
                case "Borderline SMOTE"
                    [tmp,visdata] = myBorderlineSMOTE(dataset,uniqueLabels(ii),num2Add(ii),...
                        "NumNeighbors",k, "Standardize", true);
                case "Safe-level SMOTE"
                    [tmp,visdata] = mySafeLevelSMOTE(dataset,uniqueLabels(ii),num2Add(ii),...
                        "NumNeighbors",k, "Standardize", false);
            end
            newdata = [newdata; tmp];
        end
        % updated data
        updata = [dataset; newdata];
        % pca 
        topc = table2array(updata(:,1:end-1));
        [~,score,~,~,explained,~] = pca(topc,'Centered',true,'Economy',false);
        num_comp = sum(cumsum(explained) <= 80);
        if num_comp < 2 
            num_comp = 2;
        end
        topc_red = array2table(score(:,1:num_comp));
        updata2 = addvars(topc_red, updata.label,...
            'NewVariableNames','label');
        if sqrt(size(updata2,1))*2 <= (size(updata2,2))
            [svm_res_time{si},svm_res_time_shuf{si}] = deal(nan);
            continue;
        end
        updata = updata2;
        %SVM
        SVMModel = fitcsvm(updata,'label','KernelFunction','rbf','Standardize',true,'KernelScale','auto');%
        CVSVMModel = crossval(SVMModel,'kfold', 10);
        generalizationRate = kfoldLoss(CVSVMModel);
        svm_res_time{si} = (1 - mean(generalizationRate))*100;
        % check sizes
        pred_acc = [pred_acc svm_res_time{si}];
        ftr_num = [ftr_num size(updata,2)-1];
        obs_num = [obs_num size(updata,1)];

        updata.label = updata.label(randperm(length(updata.label))); % shuffle labels
        SVMModel = fitcsvm(updata,'label','KernelFunction','rbf','Standardize',true,'KernelScale','auto');%
        CVSVMModel = crossval(SVMModel,'kfold', 10);
        generalizationRate = kfoldLoss(CVSVMModel);
        svm_res_time_shuf{si} = (1 - mean(generalizationRate))*100;
    end
    svm_pa = svm_res_time;
    svm_pa_shuf = svm_res_time_shuf;
end