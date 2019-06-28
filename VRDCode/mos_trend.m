function total_mae = mos_trend(sublist,test_idx,final_scores,y_sub,scoreMat,folder,pathVar,tvNum)

full_list = 1:43;
total_mae = struct('id',{},'abs_error',{},'offset',{},'likelihood',{});

load(strcat(folder,pathVar,'scores.mat'));
load(strcat(folder,pathVar,'zscores.mat'));

[S,E] = size(scoreMat);

MAX_ITR = 5000;     % Max number of iterations
REFRESH_RATE = 0.1;
DELTA_THR = 1e-9;

for iter=1:10
    rem_diff = [];    
    
    B = NaN(size(scoreMat));
    C = NaN(size(scoreMat));    
    
    for i=1:length(sublist)
        
        rng(length(sublist)*iter+i)
        
        %%% Sampling randomly for each subject
        vidlist = randsample(1:length(test_idx),tvNum);
        
        remlist = setdiff(test_idx,test_idx(vidlist));
        rem_idx = arrayfun(@(x)find(test_idx==x,1),remlist);        
        
        train_list = test_idx(vidlist);
        test_list = remlist;
        
        B(sublist(i),train_list) = scoreMat(sublist(i),train_list);
        C(sublist(i),test_list) = scoreMat(sublist(i),test_list);
        
        rem_diff(end+1) = nanmean((scoreMat(sublist(i),remlist))' - final_scores(rem_idx));        
    end
    
    
    %% MLE
    
    % Estimate b_s and v_s for training videos
    x_e = NaN(1,E);
    x_e(test_idx) = y_sub;
    b_s = zeros(S,1);
    v_s = sqrt(nanmean((B - repmat(x_e,S,1) - repmat(b_s,1,E)).^2 ,2));    

    likelihood = zeros(1,size(C,2));
    
    for i=1:size(B,1)
        if sum(~isnan(B(i,:))) == 0
            b_s(i) = NaN;
            v_s(i) = NaN;
        end
    end    
            
    itr = 0;
    while 1
                
        v_s_prev = v_s;
        
        % b_s
        num = B - repmat(x_e,S,1) - repmat(b_s,1,E);
        den = repmat(v_s.^2,1,E);
        dl = nansum(num./den,2);
        dl2 = nansum(-(repmat(v_s.^2,1,E).^-1),2);
        
        b_s_new = b_s - dl./dl2;
        b_s = b_s*(1 - REFRESH_RATE) + b_s_new*REFRESH_RATE;
        
        % v_s
        num = (B - repmat(x_e,S,1) - repmat(b_s,1,E)).^2 - (repmat(v_s,1,E)).^2;
        den = repmat(v_s.^3,1,E);
        dl = nansum(num./den ,2);
        num = (repmat(v_s,1,E)).^2 - 3*(B - repmat(x_e,S,1) - repmat(b_s,1,E)).^2;
        den = repmat(v_s.^4,1,E);
        dl2 = nansum(num./den ,2);
        
        v_s_new = v_s - dl./dl2;
        v_s = v_s*(1 - REFRESH_RATE) + v_s_new*REFRESH_RATE;
        
        itr = itr+1;        
        delta_v_s = nansum(abs(v_s_prev - v_s).^2)^(1/2);
        
        if delta_v_s < DELTA_THR
            break
        end
        
        if itr >= MAX_ITR
            break
        end
    end
    
    
    %% Likelihood
    for j=1:size(C,2)
        for i=1:size(C,1)
            
            if ~isnan(C(i,j))                
                    likelihood(j) = likelihood(j) + (-log(v_s(i).^2)/2) - ((C(i,j) - x_e(j) - b_s(i)).^2)/2*v_s(i).^2;                
            end
        end
    end
    likelihood_f = sum(likelihood);
    %%
    
    error = nanmean(C-(repmat(x_e,S,1) + repmat(b_s,1,E)) ,2);    
    abs_error = nanmean(abs(C-(repmat(x_e,S,1) + repmat(b_s,1,E))) ,2);     % Absolute error between predicted and raw discomfort scores   
    offset = b_s;
    abs_error(setdiff(full_list,sublist)) = [];
    error(setdiff(full_list,sublist)) = [];
    offset(setdiff(full_list,sublist)) = [];

    variance = v_s;
    variance(setdiff(full_list,sublist)) = [];
    
    
    %% MLE of b_s on test videos
    
    x_e = NaN(1,E);
    x_e(test_idx) = y_sub;
    
    b_s = zeros(S,1);
    v_s = sqrt(nanmean((C - repmat(x_e,S,1) - repmat(b_s,1,E)).^2 ,2));
    
    for i=1:size(C,1)
        if sum(~isnan(C(i,:))) == 0
            b_s(i) = NaN;
            v_s(i) = NaN;
        end
    end   
        
        
    itr = 0;
    while 1
                
        v_s_prev = v_s;
        
        % b_s
        num = C - repmat(x_e,S,1) - repmat(b_s,1,E);
        den = repmat(v_s.^2,1,E);
        dl = nansum(num./den,2);
        dl2 = nansum(-(repmat(v_s.^2,1,E).^-1),2);
        
        b_s_new = b_s - dl./dl2;
        b_s = b_s*(1 - REFRESH_RATE) + b_s_new*REFRESH_RATE;
        
        % v_s
        num = (C - repmat(x_e,S,1) - repmat(b_s,1,E)).^2 - (repmat(v_s,1,E)).^2;
        den = repmat(v_s.^3,1,E);
        dl = nansum(num./den ,2);
        num = (repmat(v_s,1,E)).^2 - 3*(C - repmat(x_e,S,1) - repmat(b_s,1,E)).^2;
        den = repmat(v_s.^4,1,E);
        dl2 = nansum(num./den ,2);
        v_s_new = v_s - dl./dl2;
        v_s = v_s*(1 - REFRESH_RATE) + v_s_new*REFRESH_RATE;
        
        itr = itr+1;
        
        delta_v_s = nansum(abs(v_s_prev - v_s).^2)^(1/2);
        
        if delta_v_s < DELTA_THR
            break
        end
        %
        if itr >= MAX_ITR
            break
        end
    end
    
    %%
    
    test_offset = b_s;      % Offset estimated on test videos
    test_offset(setdiff(full_list,sublist)) = [];
    offerr = abs(offset - test_offset);     % Offset error
    %% Offset
 
    total_mae(iter).id = iter;
    total_mae(iter).error = error';
    total_mae(iter).abs_error = abs_error';
    total_mae(iter).offset = offset';
    total_mae(iter).likelihood = likelihood_f;
    total_mae(iter).offerr = offerr';
    total_mae(iter).variance = variance';    
end

end