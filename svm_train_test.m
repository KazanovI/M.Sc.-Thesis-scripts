function [ret_svm_res_cond] = svm_train_test(all_res,cells_binned,data_act,inds,varbs,varbs_check,cts,shuff)

    m_wind = 30:35; 
    num_iter = 100;
    svm_res_cond = zeros(num_iter,2);
    activity_cts = zeros(66,3);
    for si = 1:size(data_act,2)
        dat = data_act{si};
        for ct = 1:3
            dat_z = zscore(dat(logical(inds{1,si}(:,ct)),:,:),[],[2 3]);
            dat_m = mean(mean(dat_z(:,m_wind,:),[2 3]));
            if ~isempty(dat_m)
                activity_cts(si,ct) = dat_m;
            else
                activity_cts(si,ct) = nan;
            end
        end
    end
    svm_res_cond = zeros(num_iter,2);
    for iter = 1:num_iter
        % try SVM
        for ct = 1:2
            bins = [1:3; 4:6];
            x1 = cell2mat(cellfun(@(x) x(2,:),all_res,'UniformOutput',false))'; % licks at task
            x2 = cell2mat(cellfun(@(x) x,cells_binned(cts(ct),:),'UniformOutput',false))'; % perc cells
            x3 = activity_cts(:,ct);
            Y = repelem([0 1],sum(cellfun(@(x) length(x),all_res(bins)),2));
            % mat
            X = [x1 x2 x3]; 
            X_check = X(:,varbs_check);
            nn_id = any(sum(isnan(X),2),2);
            X_check(nn_id,:)= []; % remove nan
            Y(nn_id) = [];
            % make a table with labels
            dataset = array2table(X_check,'VariableNames',varbs(2,varbs_check)); 
            label0 = repmat("Train",sum(Y==0),1);
            label1 = repmat("Test",sum(Y==1),1);
            dataset = addvars(dataset, [label0;label1],...
                'NewVariableNames','label');
            t = tabulate(dataset.label);
            uniqueLabels = string(t(:,1));
            labelCounts = cell2mat(t(:,2));
            % to which group synth data is needed
            [~,mi] = min(labelCounts);
            % how much to add
            dif_labels = abs(diff(labelCounts))-1;
            if (dif_labels == -1)
                dif_labels = 0;
            end
            if mi==1
                num2Add = [dif_labels,0];
            else
                num2Add = [0,dif_labels];
            end
            % algorithm params
            algorithm = "Safe-level SMOTE";
            k = 3;
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
                            "NumNeighbors",k, "Standardize", true);
                end
                newdata = [newdata; tmp];
            end
            % updated data
            updata = [dataset; newdata];
            if shuff 
                updata.label = updata.label(randperm(length(updata.label))); % shuffle labels
            end
    %         updata.label = updata.label(randperm(length(updata.label))); % shuffle labels
            % SVM
            SVMModel = fitcsvm(updata,'label','KernelFunction','rbf','KernelScale','auto','Standardize',true);
            CVSVMModel = crossval(SVMModel,'kfold', 10);
            generalizationRate = kfoldLoss(CVSVMModel);
            svm_res_cond(iter,ct) = (1 - mean(generalizationRate))*100;
        end
    end
    ret_svm_res_cond = svm_res_cond;
end