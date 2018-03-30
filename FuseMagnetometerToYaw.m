function [states,P] =  FuseMagnetometerToYaw(states,P,magdata)
   R_yaw = (1*pi/180)^2;
   dcm = convertQuaternion2DCM(states(1:4,1));
   [ins_roll,ins_pitch,ins_yaw]=getEulerFromDCM(dcm);
   yaw_mag  = getYawFromMag( magdata, ins_roll, ins_pitch);
   k = P(4,4)/(P(4,4)+R_yaw);
   ins_yaw = ins_yaw +k*(yaw_mag-ins_yaw);
   P = P-P*k;
end