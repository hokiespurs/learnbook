function R = rpy2dcm(roll,pitch,yaw)
% Calculate Direction Cosine Matrix Using Tait-Bryan (Z-Y-X convention)
Rx = [1 0 0; 0 cos(roll) sin(roll); 0 -sin(roll) cos(roll)];
Ry = [cos(pitch) 0 -sin(pitch); 0 1 0;sin(pitch) 0 cos(pitch)];
Rz = [cos(yaw) sin(yaw) 0; -sin(yaw) cos(yaw) 0; 0 0 1];
R = Rx*Ry*Rz;
end