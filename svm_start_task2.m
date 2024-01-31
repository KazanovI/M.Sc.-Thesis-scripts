function [ret_svm_res_cond] = svm_start_task2(all_res,perc_start,perc_task,fr_start,fr_task,varbs,varbs_check,cts,shuff)

    num_iter = 100;
    svm_res_cond = zeros(num_iter,2);

    for iter = 1:num_iter
        % try SVM
        for ct = 1:2
            x1 = cell2mat(all_res)'; % licks start/task
            x1 = x1(:);
            x2 = [cell2mat(perc_start(cts(ct),:)) cell2mat(perc_task(cts(ct),:))] ; % perc cells
            x3 = [cell2mat(cellfun(@(x) x(:,cts(ct))',fr_start,'UniformOutput',false)) ...
                cell2mat(cellfun(@(x) x(:,cts(ct))',fr_task,'UniformOutput',false))]; % fr cells
            Y = repelem([0 1],sum(cellfun(@(x) length(x),all_res),2));
            % mat
            X = [x1 x2' x3']; 
            X_check = X(:,varbs_check);
            nn_id = any(sum(isnan(X),2),2);
            X_check(nn_id,:)= 0; % remove nan
            Y(nn_id) = 0;
            % make a table with labels
            dataset = array2table(X_check,'VariableNames',varbs(2,varbs_check)); 
            label0 = repmat("Start",sum(Y==0),1);
            label1 = repmat("Task",sum(Y==1),1);
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
            % SVM
            SVMModel = fitcsvm(updata,'label','KernelFunction','rbf','KernelScale','auto','Standardize',true);
            CVSVMModel = crossval(SVMModel,'kfold', 10);
            generalizationRate = kfoldLoss(CVSVMModel);
            svm_res_cond(iter,ct) = (1 - mean(generalizationRate))*100;
        end
    end
    ret_svm_res_cond = svm_res_cond;
end