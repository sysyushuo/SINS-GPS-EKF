%% 从欧拉角中获取姿态阵--机体系到导航系
function DCM = getDCMFromEuler(roll,pitch,yaw)
DCM = [cos(pitch)*cos(yaw),sin(roll)*sin(pitch)*cos(yaw)-cos(roll)*sin(yaw),cos(roll)*sin(pitch)*cos(yaw)+sin(roll)*sin(yaw);
       cos(pitch)*sin(yaw),sin(roll)*sin(pitch)*sin(yaw)+cos(roll)*cos(yaw),cos(roll)*sin(pitch)*sin(yaw)-sin(roll)*cos(yaw);
       -sin(pitch),         sin(roll)*cos(pitch),                            cos(roll)*cos(pitch)];
end