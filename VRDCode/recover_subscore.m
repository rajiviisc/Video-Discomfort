function [final_x_e, final_b_s,final_v_s] = recover_subscore(sublist,test_idx,final_scores,folder,pathVar,flag)

full_list = 1:43;


load(strcat(folder,pathVar,'scores.mat'));
load(strcat(folder,pathVar,'zscores.mat'));

[S,E] = size(scoreMat);     % S subjects, E videos
B = NaN(size(scoreMat));


MAX_ITR = 5000;
REFRESH_RATE = 0.1;
DELTA_THR = 1e-9;

if flag == 1
    B = scoreMat;
else
    B(setdiff(full_list,sublist),:) = scoreMat(setdiff(full_list,sublist),:);
end

    
%% MLE

%% Initialisation
x_e = final_scores';        % Video discomfort score
b_s = zeros(S,1);           % Subject bias
v_s = sqrt(nanmean((B - repmat(x_e,S,1) - repmat(b_s,1,E)).^2 ,2));     % Subject unreliablility
%%

    
parfor i=1:size(B,1)
    if sum(~isnan(B(i,:))) == 0
        b_s(i) = NaN;
        v_s(i) = NaN;
    end
end

itr = 0;
while 1
    
    x_e_prev = x_e;
    
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
    
    % x_e
    num = B - repmat(x_e,S,1) - repmat(b_s,1,E);
    den = repmat(v_s.^2,1,E);
    dl = nansum(num./den);
    dl2 = nansum((-(repmat(v_s,1,E)).^2).^-1);
    x_e_new = x_e - dl./dl2;
    x_e = x_e*(1 - REFRESH_RATE) + x_e_new*REFRESH_RATE;
    
    itr = itr+1;
    delta_x_e = norm(x_e_prev - x_e);
    
    if delta_x_e < DELTA_THR
        break
    end
    
    if itr >= MAX_ITR
        break
    end
end
    

%%

final_x_e = x_e;
final_b_s = b_s';
final_v_s = v_s';


%% Confidence
expr_xe = nansum((-(repmat(final_v_s,1,E)).^2).^-1);
conf_x_e = [(final_x_e-(1.96*(sqrt(-expr_xe)).^(-1)))', (final_x_e+(1.96*(sqrt(-expr_xe)).^(-1)))'];

expr_bs = nansum(-(repmat(final_v_s'.^2,1,E).^-1),2);
conf_b_s = [final_b_s'-(1.96*(sqrt(-expr_bs)).^(-1)), final_b_s'+(1.96*(sqrt(-expr_bs)).^(-1))];

num = (repmat(final_v_s',1,E)).^2 - 3*(B-repmat(final_x_e,S,1) - repmat(final_b_s',1,E)).^2;
den = repmat(final_v_s'.^4,1,E);
expr_vs = nansum(num./den ,2);
conf_v_s = [final_v_s'-(1.96*(sqrt(-expr_bs)).^(-1)), final_v_s'+(1.96*(sqrt(-expr_vs)).^(-1))];
    

%%
end