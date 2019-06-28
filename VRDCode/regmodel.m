function [perfMat] = regmodel(features,final_scores,folder,pathVar,varargin)
    tic
    flag  = varargin{1};
    if size(varargin,2) > 1
        tvNum = varargin{2};
    end
    
    sub1_data = textscan(fopen(strcat(folder,pathVar,'Sequence_map.csv')), '%s %s %f','Delimiter',',');
    sub1_name = sub1_data{1};
    anchorList = horzcat(5:10:45,60:10:100);
    sub1_name(anchorList) = [];
    
    load(strcat(folder,pathVar,'scores.mat'));
    perfMat = struct('lcc',{},'srocc',{},'rmse',{});        % Performance metrics  
    
    lambda = logspace(-5,6,12);     % Regularization parameter
    
    win = 0.05:0.05:0.5;        % Smoothing spline parameter
    mean_vec_error_fold = struct('mse',{},'win',{});

    X = horzcat(features.tError,features.depth,features.omega,features.tVel);
    y = final_scores;        

    lin_corr = zeros(8,1);
    spe_corr = zeros(8,1);        
    rmse = zeros(8,1);
    valid_rmse = zeros(8,1);
    
    % Unique scenes in the sequence
    [unique_scenes, u1, uidx] = unique(sub1_name);
    
    for j=1:8
        discomfort_X_train = zeros(0,size(X,2));
        discomfort_X_test = zeros(0,size(X,2));
        discomfort_y_train = zeros(0,1);
        discomfort_y_test = zeros(0,1);
        test_idx = [];
        
        rng(j)
        sublist = randsample(1:43,15);
        selective_scores = calculate_mos(folder,pathVar,sublist);

        %% Test-Train Split 
        for i=1:size(unique_scenes,1)
            ids = find(i==uidx);
            for k=1:size(ids,1)
                if ismember(i,[j,mod(j+1,10),10+j,10+mod(j+1,10)])
                    test_idx(end+1) = ids(k);
                    discomfort_X_test(end+1,:) = X(ids(k),:);
                    discomfort_y_test(end+1,:) = y(ids(k),:);                
                end
            end
        end

        train_idx = setdiff(1:100, test_idx);
        %%
        
        %% Parameter Estimation
        if flag == 1
            [final_x_e, final_b_s, final_v_s] = recover_subscore(sublist,test_idx,final_scores,folder,pathVar,flag);
        else
            [final_x_e, final_b_s, final_v_s] = recover_subscore(sublist,test_idx,selective_scores,folder,pathVar,flag);
        end
        %%
        
        final_x_e = final_x_e';     % Estimated Discomfort scores
        
        discomfort_y_train(1:size(train_idx,2),:) = final_x_e(train_idx,:);
        discomfort_y_test(1:size(test_idx,2),:) = final_x_e(test_idx,:);            
        
        %% Cross-validation
        cv = cvpartition(discomfort_y_train,'KFold',5);
        for i=1:size(lambda,2)
            
            for w=1:length(win)
                tError = calc_trajError(folder,pathVar,final_scores,win(w));        % Trajectory error using given smoothing parameter                
                
                X_new = horzcat(tError,features.depth,features.omega,features.tVel);
                discomfort_X_train = X_new(train_idx,:);
                discomfort_X_test = X_new(test_idx,:);
                
                
                for l=1:cv.NumTestSets
                    
                    testindices = cv.test(l);
                    trainindices = cv.training(l);
                    
                    trainData = discomfort_X_train(trainindices,:);
                    trainScores = discomfort_y_train(trainindices);
                    testData = discomfort_X_train(testindices,:);
                    testScores = discomfort_y_train(testindices);
                    
                    b_valid = ridge(trainScores,trainData,lambda(i),0);     
                    y_valid = [ones(size(testData,1),1) testData]*b_valid;
                    
                    mse_test = rms(y_valid-testScores);
                    vec_error(l) = mse_test;
                end
                mse_win(w) = mean(vec_error);
            end
            
            [~,minw] = min(mse_win);        % Optimal smoothing parameter
            
            mean_vec_error_fold(i).mse = mse_win(minw);
            mean_vec_error_fold(i).win = win(minw);
        end

        [~,idx] = min([mean_vec_error_fold.mse]);       % Choosing best regularization hyperparameter

        %% Calculate features using chosen hyperparameters
        tError = calc_trajError(folder,pathVar,final_scores,mean_vec_error_fold(idx).win);
        X_new = horzcat(tError,features.depth,features.omega,features.tVel);           
            
        discomfort_X_train_subset = X_new(train_idx,:);
        discomfort_X_test_subset = X_new(test_idx,:);
        %%
        
        %% Ridge
        b = ridge(discomfort_y_train,discomfort_X_train_subset,lambda(idx),0);
        y_pred = [ones(size(discomfort_X_test_subset,1),1) discomfort_X_test_subset]*b;
        %%              
                        
        r = corrcoef(y_pred,discomfort_y_test);
        lin_corr(j) = r(1,2);       % LCC
        spe_corr(j) = corr(y_pred,discomfort_y_test,'Type','Spearman');     % SROCC
        rmse(j) = rms(y_pred - discomfort_y_test);        
        
        %% Mean Absolute Error and Offset for user-specific model
        if flag == 2
            total_mae = mos_trend(sublist,test_idx,selective_scores,y_pred,scoreMat,folder,pathVar,tvNum);
            

            mean_abs_error(j,:) = nanmean((reshape([total_mae.abs_error],[size(sublist,2),size(total_mae,2)]))');
            offset(j,:) = nanmean((reshape([total_mae.offset],[size(sublist,2),size(total_mae,2)]))');
            likelihood(j) = nanmean([total_mae.likelihood]);
            offerr(j,:) = nanmean((reshape([total_mae.offerr],[size(sublist,2),size(total_mae,2)]))');            
            mean_error(j,:) = nanmean((reshape([total_mae.error],[size(sublist,2),size(total_mae,2)]))');
            variance(j,:) = nanmean((reshape([total_mae.variance],[size(sublist,2),size(total_mae,2)]))');
            
        end
        %%
        
                    
                    
    end
    
    %% Performance Measures
    perfMat(1).srocc = mean(spe_corr);
    perfMat(1).lcc = mean(lin_corr);    
    perfMat(1).rmse = mean(rmse);
    perfMat(1).valid_rmse = mean(valid_rmse);

    if flag == 2
        perfMat(1).me = mean(mean_error);
        perfMat(1).mae = mean(mean_abs_error);
        perfMat(1).offset = mean(offset);
        perfMat(1).likelihood = mean(likelihood);
        perfMat(1).offerr = mean(offerr);        
        perfMat(1).subVar = mean(variance);
    end
    %%
    toc
   
end
