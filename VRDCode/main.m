%% Matlab Code for "Prediction of Discomfort due to Egomotion in Immersive Videos for
%% Virtual Reality" (ISMAR 2019)

clc
clear
close all

if ~contains(pwd,'/')
    pathVar = '\';
else
    pathVar = '/';
end

% Add current folder and all subfolders to the path
mFolder = pwd; 
addpath(genpath(mFolder));

folder = strcat(mFolder,pathVar,'Data');
full_list = 1:43;
final_scores = calculate_mos(folder,pathVar,full_list);     % MOS based on Z-scores

features = vo(folder,pathVar,final_scores);

perfMat_avgUser = regmodel(features,final_scores,folder,pathVar,1);       % Discomfort prediction for average user

%% % Discomfort prediction for specific users
perfMat_specUser = regmodel(features,final_scores,folder,pathVar,2,7);      % Specify number of training videos in the arguments   
mean_abs_error = mean(perfMat_specUser.mae);        % Mean Absolute Error
offset_error = mean(perfMat_specUser.offerr);       % Offset Error