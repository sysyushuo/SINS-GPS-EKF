%% load outdoor static data
load indoorDynamicData.mat;
% delete unuseable variables
clearvars -except IMU IMU_label MAG MAG_label EKF1 EKF1_label GPS GPS_label BARO BARO_label ;
%% find IMUStart
% data length
 [IMULength,~] = size(IMU);
 indexIMU = 1;
 [MAGLength,~] = size(MAG);
 indexMAG = 1; 
 [GPSLength,~] = size(GPS);
 indexGPS = 2;
 [BAROLength,~] = size(BARO);
 indexBARO = 1;
%% calculateIniatialtion
 % gravityNED
  gravityNED = single([0;0;9.807]);
% earthRateNED,the NED earth spin vector in rad/sec
  deg2rad = single(pi/180);
  earthRateECEF = single([0, 0, 7.2921e-005]);% ECEF坐标系下地球旋转角速度
  for i = 1:length(GPS)
   earthRateNED(i,1)  = single(cos(GPS(i,8)*deg2rad)*earthRateECEF(3));
   earthRateNED(i,2)  = single(0);
   earthRateNED(i,3)  = single(-sin(GPS(i,8)*deg2rad)*earthRateECEF(3));      
  end
  earthRateNED=earthRateNED';
%% ReadMeasurements
  %% readIMUData
    % angRate =IMU(:,3:5);
    % accel=IMU(:,6:8);
  %% readGpsData
   % ConvertGpsData(GPS_DataArrived,RefLatLongDeg,LatLongDeg,GndSpd,CourseDeg,VelD);
   deg2rad = single(pi/180);
   earthRadius = single(6378145); %地球半径
%    PosNE = single(zeros(2,1));% 2*1
%    VelNED = single(zeros(3,1));% 3*1
   GndSpd = GPS(:,12);
   CourseDeg = GPS(:,13);
   VelD= GPS(:,14);
%    LatDelta = single(zeros(1,672)); 
%    LongDelta = single(zeros(1,672));
   for i= 3:length(GPS)
    LatDelta(i)   =  GPS(i,8) - GPS(2,8);
    LongDelta(i)  =  GPS(i,9) - GPS(2,9);
    PosNE(1,i) = earthRadius * LatDelta(i)/100;% m
    PosNE(2,i) = earthRadius * cos(GPS(1,8)*deg2rad) * LongDelta (i)/100;% m
    VelNED(1,i) = GndSpd(i)*cos(CourseDeg(i)*deg2rad);% m/s
    VelNED(2,i) = GndSpd(i)*sin(CourseDeg(i)*deg2rad);% m/s
    VelNED(3,i) = VelD(i);% m/s
   end
   VelNED = VelNED';
   PosNE = PosNE';

  %% readHgtData
   Alt_BARO = BARO(:,3)-BARO(1,3);
   Alt_GPS = GPS(:,11)-GPS(2,11);

  %% readMagData
  Offset = MAG(:,6:8);
  Mag = MAG(:,3:5);
  magBias =  Offset * 0.001;
  magData = Mag* 0.001 + magBias;%机体系X/Y/Z
%% CalculateInitialStates
 %% EulAHRS
 accel = IMU(1, 6:8);
init_roll = atan2(-accel(2), -accel(3));
init_pitch = atan2(accel(1), -accel(3));

mag = MAG(indexMAG, 3:5);
hx = mag(1)*cos(init_pitch) + mag(2)*sin(init_pitch)*sin(init_roll) + mag(3)*sin(init_pitch)*cos(init_roll);
hy = mag(2)*cos(init_roll) - mag(3)*sin(init_roll);
init_yaw = atan2(-hy, hx);
if(init_yaw < 0)
    init_yaw = init_yaw + 2*pi;
end
DCM = getDCMFromEuler(init_roll, init_pitch, init_yaw);
%% states——16
 states = zeros(16,1);% 16*1
 % q0 q1 q2 q3
 q0 = 0.5*sqrt(1+DCM(1,1)+DCM(2,2)+DCM(3,3));
 q1 = 0.5*sqrt(1+DCM(1,1)-DCM(2,2)-DCM(3,3))*sign(DCM(3,2)-DCM(2,3));
 q2 = 0.5*sqrt(1-DCM(1,1)+DCM(2,2)-DCM(3,3))*sign(DCM(1,3)-DCM(3,1));
 q3 = 0.5*sqrt(1-DCM(1,1)-DCM(2,2)+DCM(3,3))*sign(DCM(2,1)-DCM(1,2));
 Quat = normalizeQuaternion([q0,q1,q2,q3]);
 % MagNED
 MagNED = DCM *magData' ;% 3*1332
 %MagNED = DCM*magData' ;
 MagNED = MagNED';
 % windVelNE;
 windVelNE = [0;0];
 % MagXYZ
 magData =  magData; 
 % PosNE
 PosNE =  PosNE;
 % VelNED
 VelNED = VelNED;
 % Alt
 Alt_GPS = Alt_GPS;
 % Delta_Angle_bias_X_Y_Z
 DeAnglebias = [0;0;0];
 % Delta_Vel_bias_X_Y_Z
 DeVelbias = [0;0;0];
%% initStates
   states(1:4,:) = Quat;
   states(5:7,:) = VelNED(1,:);
   states(8:9,:) = PosNE(1,:);
   %states(10)    = Alt_BARO(1);
   states(10)    = Alt_GPS(2);
   states(11:13,:) = DeAnglebias(1,:);
   states(14:16,:) = DeVelbias(1,:);
%    states(17:19,:) = MagNED(1,:);
%    states(20:22,:) = magData(1,:);
%% initCovariance
   P =calculateInitialData();