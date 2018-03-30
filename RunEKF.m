%% Run EKF
clear all;
clc;
%% initScript
%   initScript_outdoorStaticData;
    initScript_outdoorDynamicData;
%   initScript_indoorStaticData;
%    initScript_indoorDynamicData;
euler(:,1) = [init_roll,init_pitch,init_yaw];
velNED(:,1) = states(5:7,1);
posNED(:,1) = states(8:10,1);
%% start_run
    IMU_time = IMU(1,2); 
    angRate = IMU(1,3:5);
    accel = IMU(1,6:8);
  for i =2:IMULength 
         prevTimeIMU = IMU_time;
         IMU_time = IMU(i,2);
         dt = single((IMU_time - prevTimeIMU) / 10^6);
   % Convert IMU data to delta angles and velocities using trapezoidal integration
         prevAngRate = angRate;
         angRate=IMU(i,3:5);% 1*3
         dAng = 0.5*(angRate + prevAngRate)*dt;% 1*3
         prevAccel   = accel;
         accel = IMU(i,6:8);% 1*3
         dVel = 0.5*(accel + prevAccel)*dt; % 1*3
%% 时间更新——50HZ
  %% 状
   % UpdateStrapdownEquationsNED  
        [states,correctedDelAng,correctedDelVel,accNavMag]  = UpdateStrapdownEquationsNED(states, dAng,dVel,dt);
%         CorrectedDelAng(:,i) =a;
        Tbn = convertQuaternion2DCM(states(1:4,1));
  %% 状态协方差预测
    % CovariancePrediction  
          P  = CovariancePrediction(correctedDelAng,correctedDelVel,states,P,dt);
%  %% 量测更新 
%   %% 融合GPS的测量值VelPosNED
 % FuseVelPosNED——GPS(5HZ)
         if (indexGPS < GPSLength && GPS(indexGPS,2) > prevTimeIMU && GPS(indexGPS,2) < IMU_time)
          [states,P,innovation_GPS, varInnov_GPS]= FuseVelPosNED( states,P ,accNavMag,1,VelNED(indexGPS ,1:3),1,PosNE(indexGPS,1:2),1,Alt_GPS(indexGPS));% Fuse_Alt_GPS
          correctedDelAng   = correctedDelAng - (transpose(Tbn)*earthRateNED(:,indexGPS))'*dt;
          indexGPS = indexGPS +1;
         end

%   %% 先融合GPS测量值,再融合气压值测量值
%    FuseVelPosNED——GPS(5HZ)
%        if (indexGPS < GPSLength && GPS(indexGPS,2) > prevTimeIMU && GPS(indexGPS,2) < IMU_time)
%            [states,P,innovation_GPS, varInnov_GPS]= FuseVelPosNED( states, P  ,accNavMag,1,VelNED(indexGPS ,1:3),PosNE(indexGPS,1:2),0,0);% Fuse_GPS_VelNED_PosNE
%            indexGPS = indexGPS + 1; 
%        end
%    % FuseVelPosNED——BARO(1HZ)——only fuse Alt_BARO 
%      if (indexBARO < BAROLength && BARO (indexBARO ,2) > prevTimeIMU && BARO(indexBARO ,2) < IMU_time)
%         [states,P,innovation_GPS_BARO, varInnov_GPS_BARO]= FusePosD(states, P ,1,Alt_BARO(indexBARO));% Fuse_BARO_PosD
%         indexBARO = indexBARO + 1;
%      end
%  
%  %% 融合磁力计-10HZ
%          if(indexMAG < MAGLength && MAG(indexMAG,2) > prevTimeIMU  && MAG(indexMAG,2) < IMU_time)
%              [states,P,innovatoin_MAG,varInnov_MAG] =  FuseMagnetometer(states,P,1,magData(indexMAG,1:3),1); 
%              indexMAG = indexMAG+1;
%          end
      
 %% 融合磁力计-10HZ
%          if(indexMAG < MAGLength && MAG(indexMAG,2) > prevTimeIMU  && MAG(indexMAG,2) < IMU_time)
%              [states,P] =  FuseMagnetometerToYaw(states,P,magData(indexMAG,1:3)); 
% 
% 
%              indexMAG = indexMAG+1;
%          end
%                    
%     
        
%% 存储状态量
   Cbn = convertQuaternion2DCM(states(1:4,1));
   [roll1, pitch1, yaw1] = getEulerFromDCM(Cbn);
   [init_roll,init_pitch,init_yaw] = getEulerFromDCM( DCM);

   euler(:,i) = [roll1, pitch1, yaw1];
%    MAG_yaw(:,indexMAG)=yaw_mag ;
   velNED(:,i) = states(5:7,1);
   
   posNED(:,i) = states(8:10,1);     
   
   DeltaAngleBias(:,i) = states(11:13,1);
   
   DeltaVelBias(:,i) = states(14:16,1);    
  end
  
%% plotData
 plotData;



