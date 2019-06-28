%
clc;
clearvars;
close all;
%

if ~contains(pwd,'/')
    pathVar = '\';
else
    pathVar = '/';
end

[ParentFolderPath] = fileparts(strcat(pwd,pathVar,'full_ratings.m'));
pathName = strcat(ParentFolderPath,pathVar,'Data',pathVar,'StudyData');

A_orig = textscan(fopen(strcat(pathName,pathVar,'Subject1.csv')), '%f %f','Delimiter',',');
sub1_seq = A_orig{1};

ancList = [5:10:45, 60:10:100];
full_ratingsMat = NaN(43,110);
final_anch = [];

sub1_red = sub1_seq; sub1_red(ancList) = [];

for i=1:43    
    z = textscan(fopen(strcat(pathName,pathVar,'Subject',num2str(i),'.csv')), '%f %f','Delimiter',',');
    curr_seq = z{1};
    curr_seq(ancList) = [];
    ratings = z{2};
    diff = ratings;


    diff(ancList) = [];
    sess1 = diff(1:50); sess2 = diff(51:100);
    z_score = [(sess1-mean(sess1))./std(sess1); (sess2-mean(sess2))./std(sess2)];
    z_score = (100*(z_score+3))/6;
    
    ratings = ratings(:);

    full_ratingsMat(i,1:size(ratings,1)) = ratings';


    for j=1:size(sub1_red)
        idx = find(ismember(curr_seq,sub1_red(j)));
        zscoreMat(i,j) = z_score(idx);
        scoreMat(i,j) = diff(idx);
    end

end