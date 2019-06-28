function [rotArr,rdotArr,omegMat] = calculate_rotm(euler_mat,w_ins)
%% Calculation of camera angular velocity from Euler angles
for i=1:size(euler_mat,1)
    sx = sin(euler_mat(i,1));
    sy = sin(euler_mat(i,2));
    sz = sin(euler_mat(i,3));
    cx = cos(euler_mat(i,1));
    cy = cos(euler_mat(i,2));
    cz = cos(euler_mat(i,3));
    
    % Rotation matrix
    rotm(1,1) = cy*cz;
    rotm(1,2) = sx*sy*cz - cx*sz;
    rotm(1,3) = cx*sy*cz + sx*sz;
    rotm(2,1) = cy*sz;
    rotm(2,2) = sx*sy*sz+cx*cz;
    rotm(2,3) = cx*sy*sz-sx*cz;
    rotm(3,1) = -sy;
    rotm(3,2) = sx*cy;
    rotm(3,3) = cx*cy;
    
    % dR/dt
    rdotm(1,1) = -sz*cy*w_ins(i,3) - cz*sy*w_ins(i,2);
    rdotm(1,2) = (cz*sy*cx + sz*sx)*w_ins(i,1) + (sx*cy*cz)*w_ins(i,2) - (sx*sy*sz + cx*cz)*w_ins(i,3);
    rdotm(1,3) = (sz*cx - sx*sy*cz)*w_ins(i,1) + (cx*cy*cz)*w_ins(i,2) + (sx*cz - cx*sy*sz)*w_ins(i,3);
    rdotm(2,1) = (cz*cy)*w_ins(i,3) - (sz*sy)*w_ins(i,2);
    rdotm(2,2) = (cx*sy*sz - sx*cz)*w_ins(i,1) + (sx*cy*sz)*w_ins(i,2) + (sx*sy*cz - sz*cx)*w_ins(i,3);
    rdotm(2,3) = -(sx*sy*sz+cx*cz)*w_ins(i,1) + (cx*cy*sz)*w_ins(i,2) + (cx*sy*cz + sx*sz)*w_ins(i,3);
    rdotm(3,1) = -cy*w_ins(i,2);
    rdotm(3,2) = -sx*sy*w_ins(i,2) + cx*cy*w_ins(i,1);
    rdotm(3,3) = -cx*sy*w_ins(i,2) - sx*cy*w_ins(i,1);

    % Angular velocity vector
    omegMat(i,1) = w_ins(i,1) - sy*w_ins(i,3);
    omegMat(i,2) = cx*w_ins(i,2) + cy*sx*w_ins(i,3);
    omegMat(i,3) = cx*cy*w_ins(i,3) - sx*w_ins(i,2);
    
    rotArr(:,:,i) = rotm;
    rdotArr(:,:,i) = rdotm;
end
end