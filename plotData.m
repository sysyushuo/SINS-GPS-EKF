%% plotData
%% EKF
% roll
figure(1)
hold on;
grid on;
plot(IMU(:,2)/10^6,euler(1,:)*180/pi,'r');
plot(EKF1(:,2)/10^6,EKF1(:,3),'b');
xlabel('Time/s');ylabel('Degree/。');
legend('EKF','GCS');
% legend('EKF');
title('roll');
box on;

% pitch
figure(2)
hold on
grid on
plot(IMU(:,2)/10^6,euler(2,:)*180/pi,'r');
plot(EKF1(:,2)/10^6,EKF1(:,4),'b');
xlabel('Time/s');ylabel('Degree/。');
legend('EKF','GCS');
% legend('EKF');
title('pitch');
box on

% yaw
figure(3)
hold on
grid on
plot(IMU(:,2)/10^6,euler(3,:)*180/pi,'r');
plot(EKF1(:,2)/10^6,EKF1(:,5),'b');
xlabel('Time/s');ylabel('Degree/。');
legend('EKF','GCS');
% legend('EKF');
title('yaw');
box on

% velN
figure(4);
hold on
grid on
plot(IMU(:,2)/10^6,velNED(1,:),'r');
plot(GPS(:,2)/10^6,VelNED(:,1),'b');
xlabel('Time/s');ylabel('Vel/(m/s)');
legend('EKF','GPS');
title('VN');
box on 

% velE
figure(5)
hold on
grid on
plot(IMU(:,2)/10^6,velNED(2,:),'r');
plot(GPS(:,2)/10^6,VelNED(:,2),'b');
xlabel('Time/s');ylabel('Vel/(m/s)');
legend('EKF','GPS');
title('VE');
box on 

% velD
figure(6);
hold on
grid on
plot(IMU(:,2)/10^6,velNED(3,:),'r');
plot(GPS(:,2)/10^6,VelNED(:,3),'b');
xlabel('Time/s');ylabel('Vel/(m/s)');
legend('EKF','GPS');
title('VD');
box on 

% posN
figure(7);
hold on
grid on
plot(IMU(:,2)/10^6,posNED(1,:),'r');
plot(GPS(:,2)/10^6,PosNE(:,1),'b');
xlabel('Time/s');ylabel('Pos/m');
legend('EKF','GPS');
title('PN');
box on 

% posE
figure(8);
hold on
grid on
plot(IMU(:,2)/10^6,posNED(2,:),'r');
plot(GPS(:,2)/10^6,PosNE(:,2),'b');
xlabel('Time/s');ylabel('Pos/m');
legend('EKF','GPS');
title('PE');
box on 

% posD
figure(9);
hold on
grid on
plot(IMU(:,2)/10^6,posNED(3,:),'r');
plot(GPS(:,2)/10^6,-Alt_GPS,'b');
xlabel('Time/s');ylabel('Pos/m');
legend('EKF','GPS');
title('PD');
box on 

% deltaAngleBias
figure(10);
hold on
grid on
plot(IMU(:,2)/10^6,DeltaAngleBias(1,:),'r');
plot(IMU(:,2)/10^6,DeltaAngleBias(2,:),'b');
plot(IMU(:,2)/10^6,DeltaAngleBias(3,:),'g');
xlabel('Time/s');ylabel('Degree/。');
legend('x','y','z');
title('EKF-deltaAngleBias');
box on 

% deltaVelBias
figure(11);
hold on
grid on
plot(IMU(:,2)/10^6,DeltaVelBias(1,:),'r');
plot(IMU(:,2)/10^6,DeltaVelBias(2,:),'b');
plot(IMU(:,2)/10^6,DeltaVelBias(3,:),'g');
xlabel('Time/s');ylabel('deltaVel/(m/s)');
legend('x','y','z');
title('EKF-deltaVelBias');
box on 

%%track 
figure(12);
plot3(posNED(1,:),posNED(2,:),posNED(3,:),'r','linewidth',2);hold on;grid on
% plot3(PosNE(:,1),PosNE(:,2),Alt_GPS,'b','linewidth',2);hold on;grid on;
xlabel('m');ylabel('m');zlabel('m');
% legend('EKF','GPS');
legend('EKF');
box on;
% % magNED
% figure(12);
% hold on
% grid on
% plot(IMU(:,2)/10^6,magNED_n(1,:),'r');
% plot(IMU(:,2)/10^6,magNED_n(2,:),'b');
% plot(IMU(:,2)/10^6,magNED_n(3,:),'g');
% xlabel('Time/s');ylabel('milliGaus');
% legend('N','E','D');
% title('EKFMagFluxMea！！NED');
% box on

% % magXYZ
% figure(13);
% hold on
% grid on
% plot(IMU(:,2)/10^6,magData_b(1,:),'r');
% plot(IMU(:,2)/10^6,magData_b(2,:),'b');
% plot(IMU(:,2)/10^6,magData_b(3,:),'g');
% xlabel('Time/s');ylabel('milliGaus');
% legend('x','y','z');
% title('EKFMagFluxMea！！XYZ');
% box on

% % magNED-without processing
% figure(14);
% hold on
% grid on
% plot(MAG(:,2)/10^6,MagNED(:,1),'r');
% plot(MAG(:,2)/10^6,MagNED(:,2),'b');
% plot(MAG(:,2)/10^6,MagNED(:,3),'g');
% xlabel('Time/s');ylabel('milliGaus');
% legend('N','E','D');
% title('MagFluxMea！！NED');
% box on

% % magXYZ-without processing
% figure(15);
% hold on
% grid on
% plot(MAG(:,2)/10^6,magData(:,1),'r');
% plot(MAG(:,2)/10^6,magData(:,2),'b');
% plot(MAG(:,2)/10^6,magData(:,3),'g');
% legend('x','y','z');
% title('MagFluxMea！！XYZ');
% box on
% posNED

% figure(17)
% plot(IMU(:,2)/10^6,MAG_yaw*180/pi);



%% SINS
% figure;
% hold on;
% grid on;
% subplot(3,1,1);plot(IMU(:,2)/10^6,euler(1,:)*180/pi,'r');xlabel('Time/s');ylabel('Degree/。');legend('roll');title('roll');box on;
% subplot(3,1,2);plot(IMU(:,2)/10^6,euler(2,:)*180/pi,'g');xlabel('Time/s');ylabel('Degree/。');legend('pitch');title('pitch');box on;
% subplot(3,1,3);plot(IMU(:,2)/10^6,euler(3,:)*180/pi,'b');xlabel('Time/s');ylabel('Degree/。');legend('yaw');title('yaw');box on;
% subplot(3,1,1);plot(IMU(:,2)/10^6,velNED(1,:),'r');xlabel('Time/s');ylabel('Vel/(m/s)');legend('Vn');title('velcity！North');box on; 
% subplot(3,1,2 );plot(IMU(:,2)/10^6,velNED(2,:),'g');xlabel('Time/s');ylabel('Vel/(m/s)');legend('Ve');title('velcity！East');box on; 
% subplot(3,1,3);plot(IMU(:,2)/10^6,velNED(3,:),'b');xlabel('Time/s');ylabel('Vel/(m/s)');legend('Vd');title('velcity！Down');box on; 
% subplot(3,1,1);plot(IMU(:,2)/10^6,posNED(1,:),'r');xlabel('Time/s');ylabel('Pos/m');legend('Pn');title('position！North');box on;
% subplot(3,1,2);plot(IMU(:,2)/10^6,posNED(2,:),'g');xlabel('Time/s');ylabel('Pos/m');legend('Pe');title('position！East');box on;
% subplot(3,1,3);plot(IMU(:,2)/10^6,posNED(3,:),'b');xlabel('Time/s');ylabel('Pos/m');legend('Pd');title('position！Down');box on;


%% single Sample 
% figure(1);
% hold on;
% grid on;
% subplot(3,1,1);plot(IMU(:,2)/10^6,CorrectedDelAng(1,:),'r');title('Single sample');legend('theta-x');xlabel('Time/s');ylabel('rad/s');
% subplot(3,1,2);plot(IMU(:,2)/10^6,CorrectedDelAng(2,:),'g');title('Single sample');legend('theta-y');xlabel('Time/s');ylabel('rad/s');
% subplot(3,1,3);plot(IMU(:,2)/10^6,CorrectedDelAng(3,:),'b');title('Single sample');legend('theta-z');xlabel('Time/s');ylabel('rad/s');

%% earth rotation
% figure(1);
% hold on;
% grid on;
% subplot(3,1,1);plot(GPS(:,2)/10^6,earthRateNED(1,:),'r');title('earth rotation');legend('earth-E');xlabel('Time/s');ylabel('degree/s');   
% subplot(3,1,2);plot(GPS(:,2)/10^6,earthRateNED(2,:),'g');title('earth rotation');legend('earth-N');xlabel('Time/s');ylabel('degree/s');   
% subplot(3,1,3);plot(GPS(:,2)/10^6,earthRateNED(3,:),'b');title('earth rotation');legend('earth-D');xlabel('Time/s');ylabel('degree/s');   

%% rotational correction
% figure(1);
% hold on;
% grid on;
% subplot(3,1,1);plot(IMU(:,2)/10^6,CorrectedDelAng(1,:),'r');title('rotational correction');legend('delta-Vx'); xlabel('Time/s');ylabel('m/s');
% subplot(3,1,2);plot(IMU(:,2)/10^6,CorrectedDelAng(2,:),'g');title('rotational correction');legend('delta-Vy'); xlabel('Time/s');ylabel('m/s');
% subplot(3,1,3);plot(IMU(:,2)/10^6,CorrectedDelAng(3,:),'b');title('rotational correction');legend('delta-Vz'); xlabel('Time/s');ylabel('m/s');

%% skulling correction
% figure(1);
% hold on;
% grid on;
% subplot(3,1,1);plot(IMU(:,2)/10^6,CorrectedDelAng(1,:),'r');title('skulling correction');legend('delta-Vx'); xlabel('Time/s');ylabel('m/s');
% subplot(3,1,2);plot(IMU(:,2)/10^6,CorrectedDelAng(2,:),'g');title('skulling correction');legend('delta-Vy'); xlabel('Time/s');ylabel('m/s');
% subplot(3,1,3);plot(IMU(:,2)/10^6,CorrectedDelAng(3,:),'b');title('skulling correction');legend('delta-Vz'); xlabel('Time/s');ylabel('m/s');

%% GPS_VelNED、GPSPosNED
% figure(1);
% subplot(3,1,1);plot(GPS(:,2)/10^6,VelNED(:,1),'r');title('GPS-VE');legend('GPSVelNED');xlabel('Time/s');ylabel('m/s');
% subplot(3,1,2);plot(GPS(:,2)/10^6,VelNED(:,2),'g');title('GPS-VN');legend('GPSVelNED');xlabel('Time/s');ylabel('m/s'); 
% subplot(3,1,3);plot(GPS(:,2)/10^6,VelNED(:,3),'b');title('GPS-VD');legend('GPSVelNED');xlabel('Time/s');ylabel('m/s'); 
% figure(2);
% subplot(3,1,1);plot(GPS(:,2)/10^6,PosNE(:,1),'r');title('GPS-PE');legend('GPSPosNED');xlabel('Time/s');ylabel('m');
% subplot(3,1,2);plot(GPS(:,2)/10^6,PosNE(:,2),'g');title('GPS-PN');legend('GPSPosNED');xlabel('Time/s');ylabel('m'); 
% subplot(3,1,3);plot(GPS(:,2)/10^6,Alt_GPS,'b');title('GPS-PD');legend('GPSPosNED');xlabel('Time/s');ylabel('m'); 