function tError = calc_trajError(pathName,pathVar,final_scores,wId)

A_orig = textscan(fopen(strcat(pathName,pathVar,'StudyData',pathVar,'Subject1.csv')), '%f %f','Delimiter',',');
videoId = A_orig{1};
sub1_scores = A_orig{2};
fclose('all');

frames = (1:300)';
anchorList = horzcat(5:10:45,60:10:100);
finalVidList = videoId;
finalVidList(anchorList) = [];

meanFilterFunction = @(theBlockStructure) mean(theBlockStructure.data(:));

for vid=1:length(finalVidList) 
    
    Pose = table2array(readtable(char(strcat(pathName,pathVar,'Camera Pose',pathVar,'Video',num2str(finalVidList(vid)),'.txt'))));
    T = Pose(:,1:3);
    R = Pose(:,4:6);
    
    if ismember(finalVidList(vid),[14,49,57,73,76,97])
        T = T.*100;
    end
    
    f{1} = csaps(frames,T(:,1),wId); f{2} = csaps(frames,T(:,2),wId); f{3} = csaps(frames,T(:,3),wId);
    f{4} = csaps(frames,R(:,1),wId); f{5} = csaps(frames,R(:,2),wId); f{6} = csaps(frames,R(:,3),wId);
    for smoothVar = 1:size(T,2)
        T_smooth(:,smoothVar) = ppval(f{1,smoothVar},frames);
        R_smooth(:,smoothVar) = ppval(f{1,smoothVar+3},frames);
    end
    
    t_error = abs(T-T_smooth);
    r_error = abs(R-R_smooth);    
    
    blockMeans = blockproc([t_error r_error], [30,1], meanFilterFunction);
    tError(vid,:) = reshape(blockMeans', [1,size(blockMeans,1)*size(blockMeans,2)]);
end

end