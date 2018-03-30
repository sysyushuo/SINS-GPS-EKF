%% convert quaternion to dcm
% body to reference
function DCM = convertQuaternion2DCM(quat)

DCM = zeros(3,3);
q0 = quat(1);
q1 = quat(2);
q2 = quat(3);
q3 = quat(4);
DCM(1,1) = q1^2 + q0^2 - q3^2 - q2^2;
DCM(1,2) = 2*(q1*q2 - q0*q3);
DCM(1,3) = 2*(q1*q3 + q0*q2);
DCM(2,1) = 2*(q1*q2 + q0*q3);
DCM(2,2) = q2^2 - q3^2 + q0^2 - q1^2;
DCM(2,3) = 2*(q2*q3 - q0*q1);
DCM(3,1) = 2*(q1*q3 - q0*q2);
DCM(3,2) = 2*(q2*q3 + q0*q1);
DCM(3,3) = q3^2 - q2^2 - q1^2 + q0^2;

end