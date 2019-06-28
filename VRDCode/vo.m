function features = vo(pathName,pathVar,final_scores)
A_orig = textscan(fopen(strcat(pathName,pathVar,'StudyData',pathVar,'Subject1.csv')), '%f %f','Delimiter',',');
videoId = A_orig{1};
fclose('all');

avg_vel = zeros(100,30);
frames = (1:300)';

% Positions of anchor videos in each sequence
anchorList = horzcat(5:10:45,60:10:100);
finalVidList = videoId;
finalVidList(anchorList) = [];

wId = 0.05;      % Smoothing spline parameter
avg_omega = zeros(100,30);

load(strcat(pathName,pathVar,'avg_depth.mat'));

meanFilterFunction = @(theBlockStructure) mean(theBlockStructure.data(:));


for vid=1:length(finalVidList)    
    
    Pose = table2array(readtable(char(strcat(pathName,pathVar,'Camera Pose',pathVar,'Video',num2str(finalVidList(vid)),'.txt'))));
    
    % T and R are the translation and euler angles matrices respectively
    T = Pose(:,1:3);
    R = Pose(:,4:6);
    
    % Scale of translation values readjusted to match the actual scale used
    % in the scene
    if ismember(finalVidList(vid),[14,49,57,73,76,97])
        T = T.*100;
    end
            
    X_rot = R(:,1); Y_rot = R(:,2); Z_rot = R(:,3); X_tr = T(:,1); Y_tr = T(:,2); Z_tr = T(:,3);
    
    % ppform of cubic smoothing spline 
    f{1} = csaps(frames,X_tr,wId); f{2} = csaps(frames,Y_tr,wId); f{3} = csaps(frames,Z_tr,wId);
    f{4} = csaps(frames,X_rot,wId); f{5} = csaps(frames,Y_rot,wId); f{6} = csaps(frames,Z_rot,wId);
    
    for smoothVar = 1:size(T,2)
        T_smooth(:,smoothVar) = ppval(f{1,smoothVar},frames);
        R_smooth(:,smoothVar) = ppval(f{1,smoothVar+3},frames);
    end
    
    % Instantaneous linear and angular velocities    
    v_ins = [gradient(X_tr) gradient(Y_tr) gradient(Z_tr)];
    w_ins = [gradient(X_rot) gradient(Y_rot) gradient(Z_rot)];    
    
    
    % Angular velocity of camera from euler angles   
    euler_ang = deg2rad(R(:,1:3));
    [~,~,omegaMat] = calculate_rotm(euler_ang,w_ins);

    aw = blockproc(abs(omegaMat), [30,1], meanFilterFunction);
    avg_omega(vid,:) = reshape(aw',[1 30]);

    % Compute absolute trajectory error
    t_error = abs(T-T_smooth);
    r_error = abs(R-R_smooth);
    
    % Avg trajectory error     
    avg_error = blockproc([t_error r_error], [30,1], meanFilterFunction);    
    
    % Avg absolute linear velocity
    avg_trans = blockproc(abs(v_ins), [30,1], meanFilterFunction);    
   
    trajectory_error(vid,:) = reshape(avg_error', [1,size(avg_error,1)*size(avg_error,2)]);    
    rotError(vid,:) = reshape((avg_error(:,4:6))', [1,size(avg_error(:,4:6),1)*size(avg_error(:,4:6),2)]);
    transError(vid,:) = reshape((avg_error(:,1:3))', [1,size(avg_error(:,1:3),1)*size(avg_error(:,1:3),2)]);
    
    % Average translational velocity at different scales
    avg_vel(vid,:) = reshape(avg_trans', [1,30]);

end

features.tError = trajectory_error;                                  
features.transError = transError;
features.rotError = rotError;
features.tVel = avg_vel;                                            
features.depth = avg_depth;
features.omega = avg_omega;                                         
fclose('all');