%#codegen
function [...
    states, ... % state vector after fusion of measurements
    P, ... % state covariance matrix after fusion of corrections
    innovation, ... % Magnetometer innovations - milligauss
    varInnov] ... % Innovation variances - millgauss^2
    = FuseMagnetometer( ...
    states, ... % predicted states
    P, ... % predicted covariance
    FuseData, ... % Boolean to initiate fusion of magnetometer measurements
    MagData, ... % Magnetic flux measurements (milliGauss)
    useCompass) % Boolean to indicate that a Magnetometer is being used

persistent q0 q1 q2 q3 magN magE magD;
if isempty(q0)
    q0 = single(1);
end
if isempty(q1)
    q1 = single(0);
end
if isempty(q2)
    q2 = single(0);
end
if isempty(q3)
    q3 = single(0);
end
if isempty(magN)
    magN = single(0);
end
if isempty(magE)
    magE = single(0);
end
if isempty(magD)
    magD = single(0);
end
persistent obsIndex
if isempty(obsIndex)
    obsIndex = 1;
end
persistent SH_MAG
if isempty(SH_MAG)
    SH_MAG = single(zeros(1,9));
end
persistent Tnb
if isempty(Tnb)
    Tnb = single(eye(3,3));
end
persistent MagPred
if isempty(MagPred)
    MagPred = single(zeros(3,1));
end

% declare arrays used to calculate intermediate result in covariance
% and state update
Kfusion = single(zeros(22,1));
H_MAG = single(zeros(1,22));
KH    = single(zeros(22,22));
KHP = single(zeros(22,22));
% define magnetometer measurement error variance (milligauss).
R_MAG = single(50^2);
% declare innovation variances
varInnov = single(zeros(3,1));
% declare innovation vector
innovation = single(zeros(1,3));
% Perform sequential fusion of Magnetometer measurements.
% This assumes that the errors in the dfferent componenets are
% uncorrelated which is not true, however in the absence of covariance
% data fit is the only assumption we can make
% so we might as well take advantage of the computational efficiencies
% associated with sequential fusion
if (useCompass && FuseData) || (useCompass && (obsIndex == 2 || obsIndex == 3))
    % Sequential fusion of XYZ components to spread processing load across
    % three prediction time steps.
    % Calculate observation jacobians and Kalman gains
    if FuseData
        % Copy required states to local variable names
        q0 =states(1);
        q1 =states(2);
        q2 =states(3);
        q3 =states(4);
        magN =states(17);
        magE =states(18);
        magD =states(19);
        % rotate predicted earth components into body axes and calculate
        % predicted measurments
        Tnb = single([q0^2 + q1^2 - q2^2 - q3^2, 2*(q1*q2 + q0*q3) , 2*(q1*q3-q0*q2) ;...
            2*(q1*q2 - q0*q3), q0^2 - q1^2 + q2^2 - q3^2, 2*(q2*q3 + q0*q1);...
            2*(q1*q3 + q0*q2) , 2*(q2*q3 - q0*q1) , q0^2 - q1^2 - q2^2 + q3^2]);
        MagPred = Tnb*states(17:19) +states(20:22);
        % Calculate observation jacobian common terms
        SH_MAG(1) = 2*magD*q3 + 2*magE*q2 + 2*magN*q1;
        SH_MAG(2) = 2*magD*q0 - 2*magE*q1 + 2*magN*q2;
        SH_MAG(3) = 2*magD*q1 + 2*magE*q0 - 2*magN*q3;
        SH_MAG(4) = q3^2;
        SH_MAG(5) = q2^2;
        SH_MAG(6) = q1^2;
        SH_MAG(7) = q0^2;
        SH_MAG(8) = 2*magN*q0;
        SH_MAG(9) = 2*magE*q3;
        % calculate observation jacobian
        H_MAG(1) = SH_MAG(8) + SH_MAG(9) - 2*magD*q2;
        H_MAG(2) = SH_MAG(1);
        H_MAG(3) = -SH_MAG(2);
        H_MAG(4) = SH_MAG(3);
        H_MAG(17) = SH_MAG(6) - SH_MAG(5) - SH_MAG(4) + SH_MAG(7);
        H_MAG(18) = 2*q0*q3 + 2*q1*q2;
        H_MAG(19) = 2*q1*q3 - 2*q0*q2;
        H_MAG(20) = 1;
        % calculate Kalman gain
        SK_MX = zeros(5,1);
        SK_MX(1) = 1/(P(20,20) + R_MAG + P(2,20)*SH_MAG(1) - P(3,20)*SH_MAG(2) + P(4,20)*SH_MAG(3) - P(17,20)*(SH_MAG(4) + SH_MAG(5) - SH_MAG(6) - SH_MAG(7)) + (2*q0*q3 + 2*q1*q2)*(P(20,18) + P(2,18)*SH_MAG(1) - P(3,18)*SH_MAG(2) + P(4,18)*SH_MAG(3) - P(17,18)*(SH_MAG(4) + SH_MAG(5) - SH_MAG(6) - SH_MAG(7)) + P(18,18)*(2*q0*q3 + 2*q1*q2) - P(19,18)*(2*q0*q2 - 2*q1*q3) + P(1,18)*(SH_MAG(8) + SH_MAG(9) - 2*magD*q2)) - (2*q0*q2 - 2*q1*q3)*(P(20,19) + P(2,19)*SH_MAG(1) - P(3,19)*SH_MAG(2) + P(4,19)*SH_MAG(3) - P(17,19)*(SH_MAG(4) + SH_MAG(5) - SH_MAG(6) - SH_MAG(7)) + P(18,19)*(2*q0*q3 + 2*q1*q2) - P(19,19)*(2*q0*q2 - 2*q1*q3) + P(1,19)*(SH_MAG(8) + SH_MAG(9) - 2*magD*q2)) + (SH_MAG(8) + SH_MAG(9) - 2*magD*q2)*(P(20,1) + P(2,1)*SH_MAG(1) - P(3,1)*SH_MAG(2) + P(4,1)*SH_MAG(3) - P(17,1)*(SH_MAG(4) + SH_MAG(5) - SH_MAG(6) - SH_MAG(7)) + P(18,1)*(2*q0*q3 + 2*q1*q2) - P(19,1)*(2*q0*q2 - 2*q1*q3) + P(1,1)*(SH_MAG(8) + SH_MAG(9) - 2*magD*q2)) + P(18,20)*(2*q0*q3 + 2*q1*q2) - P(19,20)*(2*q0*q2 - 2*q1*q3) + SH_MAG(1)*(P(20,2) + P(2,2)*SH_MAG(1) - P(3,2)*SH_MAG(2) + P(4,2)*SH_MAG(3) - P(17,2)*(SH_MAG(4) + SH_MAG(5) - SH_MAG(6) - SH_MAG(7)) + P(18,2)*(2*q0*q3 + 2*q1*q2) - P(19,2)*(2*q0*q2 - 2*q1*q3) + P(1,2)*(SH_MAG(8) + SH_MAG(9) - 2*magD*q2)) - SH_MAG(2)*(P(20,3) + P(2,3)*SH_MAG(1) - P(3,3)*SH_MAG(2) + P(4,3)*SH_MAG(3) - P(17,3)*(SH_MAG(4) + SH_MAG(5) - SH_MAG(6) - SH_MAG(7)) + P(18,3)*(2*q0*q3 + 2*q1*q2) - P(19,3)*(2*q0*q2 - 2*q1*q3) + P(1,3)*(SH_MAG(8) + SH_MAG(9) - 2*magD*q2)) + SH_MAG(3)*(P(20,4) + P(2,4)*SH_MAG(1) - P(3,4)*SH_MAG(2) + P(4,4)*SH_MAG(3) - P(17,4)*(SH_MAG(4) + SH_MAG(5) - SH_MAG(6) - SH_MAG(7)) + P(18,4)*(2*q0*q3 + 2*q1*q2) - P(19,4)*(2*q0*q2 - 2*q1*q3) + P(1,4)*(SH_MAG(8) + SH_MAG(9) - 2*magD*q2)) - (SH_MAG(4) + SH_MAG(5) - SH_MAG(6) - SH_MAG(7))*(P(20,17) + P(2,17)*SH_MAG(1) - P(3,17)*SH_MAG(2) + P(4,17)*SH_MAG(3) - P(17,17)*(SH_MAG(4) + SH_MAG(5) - SH_MAG(6) - SH_MAG(7)) + P(18,17)*(2*q0*q3 + 2*q1*q2) - P(19,17)*(2*q0*q2 - 2*q1*q3) + P(1,17)*(SH_MAG(8) + SH_MAG(9) - 2*magD*q2)) + P(1,20)*(SH_MAG(8) + SH_MAG(9) - 2*magD*q2));
        SK_MX(2) = SH_MAG(4) + SH_MAG(5) - SH_MAG(6) - SH_MAG(7);
        SK_MX(3) = SH_MAG(8) + SH_MAG(9) - 2*magD*q2;
        SK_MX(4) = 2*q0*q2 - 2*q1*q3;
        SK_MX(5) = 2*q0*q3 + 2*q1*q2;
        

        Kfusion(1) = SK_MX(1)*(P(1,20) + P(1,2)*SH_MAG(1) - P(1,3)*SH_MAG(2) + P(1,4)*SH_MAG(3) + P(1,1)*SK_MX(3) - P(1,17)*SK_MX(2) + P(1,18)*SK_MX(5) - P(1,19)*SK_MX(4));
        Kfusion(2) = SK_MX(1)*(P(2,20) + P(2,2)*SH_MAG(1) - P(2,3)*SH_MAG(2) + P(2,4)*SH_MAG(3) + P(2,1)*SK_MX(3) - P(2,17)*SK_MX(2) + P(2,18)*SK_MX(5) - P(2,19)*SK_MX(4));
        Kfusion(3) = SK_MX(1)*(P(3,20) + P(3,2)*SH_MAG(1) - P(3,3)*SH_MAG(2) + P(3,4)*SH_MAG(3) + P(3,1)*SK_MX(3) - P(3,17)*SK_MX(2) + P(3,18)*SK_MX(5) - P(3,19)*SK_MX(4));
        Kfusion(4) = SK_MX(1)*(P(4,20) + P(4,2)*SH_MAG(1) - P(4,3)*SH_MAG(2) + P(4,4)*SH_MAG(3) + P(4,1)*SK_MX(3) - P(4,17)*SK_MX(2) + P(4,18)*SK_MX(5) - P(4,19)*SK_MX(4));
        Kfusion(5) = SK_MX(1)*(P(5,20) + P(5,2)*SH_MAG(1) - P(5,3)*SH_MAG(2) + P(5,4)*SH_MAG(3) + P(5,1)*SK_MX(3) - P(5,17)*SK_MX(2) + P(5,18)*SK_MX(5) - P(5,19)*SK_MX(4));
        Kfusion(6) = SK_MX(1)*(P(6,20) + P(6,2)*SH_MAG(1) - P(6,3)*SH_MAG(2) + P(6,4)*SH_MAG(3) + P(6,1)*SK_MX(3) - P(6,17)*SK_MX(2) + P(6,18)*SK_MX(5) - P(6,19)*SK_MX(4));
        Kfusion(7) = SK_MX(1)*(P(7,20) + P(7,2)*SH_MAG(1) - P(7,3)*SH_MAG(2) + P(7,4)*SH_MAG(3) + P(7,1)*SK_MX(3) - P(7,17)*SK_MX(2) + P(7,18)*SK_MX(5) - P(7,19)*SK_MX(4));
        Kfusion(8) = SK_MX(1)*(P(8,20) + P(8,2)*SH_MAG(1) - P(8,3)*SH_MAG(2) + P(8,4)*SH_MAG(3) + P(8,1)*SK_MX(3) - P(8,17)*SK_MX(2) + P(8,18)*SK_MX(5) - P(8,19)*SK_MX(4));
        Kfusion(9) = SK_MX(1)*(P(9,20) + P(9,2)*SH_MAG(1) - P(9,3)*SH_MAG(2) + P(9,4)*SH_MAG(3) + P(9,1)*SK_MX(3) - P(9,17)*SK_MX(2) + P(9,18)*SK_MX(5) - P(9,19)*SK_MX(4));
        Kfusion(10) = SK_MX(1)*(P(10,20) + P(10,2)*SH_MAG(1) - P(10,3)*SH_MAG(2) + P(10,4)*SH_MAG(3) + P(10,1)*SK_MX(3) - P(10,17)*SK_MX(2) + P(10,18)*SK_MX(5) - P(10,19)*SK_MX(4));
        Kfusion(11) = SK_MX(1)*(P(11,20) + P(11,2)*SH_MAG(1) - P(11,3)*SH_MAG(2) + P(11,4)*SH_MAG(3) + P(11,1)*SK_MX(3) - P(11,17)*SK_MX(2) + P(11,18)*SK_MX(5) - P(11,19)*SK_MX(4));
        Kfusion(12) = SK_MX(1)*(P(12,20) + P(12,2)*SH_MAG(1) - P(12,3)*SH_MAG(2) + P(12,4)*SH_MAG(3) + P(12,1)*SK_MX(3) - P(12,17)*SK_MX(2) + P(12,18)*SK_MX(5) - P(12,19)*SK_MX(4));
        Kfusion(13) = SK_MX(1)*(P(13,20) + P(13,2)*SH_MAG(1) - P(13,3)*SH_MAG(2) + P(13,4)*SH_MAG(3) + P(13,1)*SK_MX(3) - P(13,17)*SK_MX(2) + P(13,18)*SK_MX(5) - P(13,19)*SK_MX(4));
        Kfusion(14) = SK_MX(1)*(P(14,20) + P(14,2)*SH_MAG(1) - P(14,3)*SH_MAG(2) + P(14,4)*SH_MAG(3) + P(14,1)*SK_MX(3) - P(14,17)*SK_MX(2) + P(14,18)*SK_MX(5) - P(14,19)*SK_MX(4));
        Kfusion(15) = SK_MX(1)*(P(15,20) + P(15,2)*SH_MAG(1) - P(15,3)*SH_MAG(2) + P(15,4)*SH_MAG(3) + P(15,1)*SK_MX(3) - P(15,17)*SK_MX(2) + P(15,18)*SK_MX(5) - P(15,19)*SK_MX(4));
        Kfusion(16) = SK_MX(1)*(P(16,20) + P(16,2)*SH_MAG(1) - P(16,3)*SH_MAG(2) + P(16,4)*SH_MAG(3) + P(16,1)*SK_MX(3) - P(16,17)*SK_MX(2) + P(16,18)*SK_MX(5) - P(16,19)*SK_MX(4));
        Kfusion(17) = SK_MX(1)*(P(17,20) + P(17,2)*SH_MAG(1) - P(17,3)*SH_MAG(2) + P(17,4)*SH_MAG(3) + P(17,1)*SK_MX(3) - P(17,17)*SK_MX(2) + P(17,18)*SK_MX(5) - P(17,19)*SK_MX(4));
        Kfusion(18) = SK_MX(1)*(P(18,20) + P(18,2)*SH_MAG(1) - P(18,3)*SH_MAG(2) + P(18,4)*SH_MAG(3) + P(18,1)*SK_MX(3) - P(18,17)*SK_MX(2) + P(18,18)*SK_MX(5) - P(18,19)*SK_MX(4));
        Kfusion(19) = SK_MX(1)*(P(19,20) + P(19,2)*SH_MAG(1) - P(19,3)*SH_MAG(2) + P(19,4)*SH_MAG(3) + P(19,1)*SK_MX(3) - P(19,17)*SK_MX(2) + P(19,18)*SK_MX(5) - P(19,19)*SK_MX(4));
        Kfusion(20) = SK_MX(1)*(P(20,20) + P(20,2)*SH_MAG(1) - P(20,3)*SH_MAG(2) + P(20,4)*SH_MAG(3) + P(20,1)*SK_MX(3) - P(20,17)*SK_MX(2) + P(20,18)*SK_MX(5) - P(20,19)*SK_MX(4));
        Kfusion(21) = SK_MX(1)*(P(21,20) + P(21,2)*SH_MAG(1) - P(21,3)*SH_MAG(2) + P(21,4)*SH_MAG(3) + P(21,1)*SK_MX(3) - P(21,17)*SK_MX(2) + P(21,18)*SK_MX(5) - P(21,19)*SK_MX(4));
        Kfusion(22) = SK_MX(1)*(P(22,20) + P(22,2)*SH_MAG(1) - P(22,3)*SH_MAG(2) + P(22,4)*SH_MAG(3) + P(22,1)*SK_MX(3) - P(22,17)*SK_MX(2) + P(22,18)*SK_MX(5) - P(22,19)*SK_MX(4));
        varInnov(1) = 1/SK_MX(1);
        % reset the observation index to 1 (we start by fusing the X
        % measurement)
        obsIndex = 1;
    elseif obsIndex == 2 % we are now fusing the Y measurement
        H_MAG(1) = SH_MAG(3);
        H_MAG(2) = SH_MAG(2);
        H_MAG(3) = SH_MAG(1);
        H_MAG(4) = 2*magD*q2 - SH_MAG(9) - SH_MAG(8);
        H_MAG(17) = 2*q1*q2 - 2*q0*q3;
        H_MAG(18) = SH_MAG(5) - SH_MAG(4) - SH_MAG(6) + SH_MAG(7);
        H_MAG(19) = 2*q0*q1 + 2*q2*q3;
        H_MAG(21) = 1;
        % calculate the Kalman gain
        SK_MY = zeros(5,1);
        SK_MY(1) = 1/(P(21,21) + R_MAG + P(1,21)*SH_MAG(3) + P(2,21)*SH_MAG(2) + P(3,21)*SH_MAG(1) - P(18,21)*(SH_MAG(4) - SH_MAG(5) + SH_MAG(6) - SH_MAG(7)) - (2*q0*q3 - 2*q1*q2)*(P(21,17) + P(1,17)*SH_MAG(3) + P(2,17)*SH_MAG(2) + P(3,17)*SH_MAG(1) - P(18,17)*(SH_MAG(4) - SH_MAG(5) + SH_MAG(6) - SH_MAG(7)) - P(17,17)*(2*q0*q3 - 2*q1*q2) + P(19,17)*(2*q0*q1 + 2*q2*q3) - P(4,17)*(SH_MAG(8) + SH_MAG(9) - 2*magD*q2)) + (2*q0*q1 + 2*q2*q3)*(P(21,19) + P(1,19)*SH_MAG(3) + P(2,19)*SH_MAG(2) + P(3,19)*SH_MAG(1) - P(18,19)*(SH_MAG(4) - SH_MAG(5) + SH_MAG(6) - SH_MAG(7)) - P(17,19)*(2*q0*q3 - 2*q1*q2) + P(19,19)*(2*q0*q1 + 2*q2*q3) - P(4,19)*(SH_MAG(8) + SH_MAG(9) - 2*magD*q2)) - (SH_MAG(8) + SH_MAG(9) - 2*magD*q2)*(P(21,4) + P(1,4)*SH_MAG(3) + P(2,4)*SH_MAG(2) + P(3,4)*SH_MAG(1) - P(18,4)*(SH_MAG(4) - SH_MAG(5) + SH_MAG(6) - SH_MAG(7)) - P(17,4)*(2*q0*q3 - 2*q1*q2) + P(19,4)*(2*q0*q1 + 2*q2*q3) - P(4,4)*(SH_MAG(8) + SH_MAG(9) - 2*magD*q2)) - P(17,21)*(2*q0*q3 - 2*q1*q2) + P(19,21)*(2*q0*q1 + 2*q2*q3) + SH_MAG(3)*(P(21,1) + P(1,1)*SH_MAG(3) + P(2,1)*SH_MAG(2) + P(3,1)*SH_MAG(1) - P(18,1)*(SH_MAG(4) - SH_MAG(5) + SH_MAG(6) - SH_MAG(7)) - P(17,1)*(2*q0*q3 - 2*q1*q2) + P(19,1)*(2*q0*q1 + 2*q2*q3) - P(4,1)*(SH_MAG(8) + SH_MAG(9) - 2*magD*q2)) + SH_MAG(2)*(P(21,2) + P(1,2)*SH_MAG(3) + P(2,2)*SH_MAG(2) + P(3,2)*SH_MAG(1) - P(18,2)*(SH_MAG(4) - SH_MAG(5) + SH_MAG(6) - SH_MAG(7)) - P(17,2)*(2*q0*q3 - 2*q1*q2) + P(19,2)*(2*q0*q1 + 2*q2*q3) - P(4,2)*(SH_MAG(8) + SH_MAG(9) - 2*magD*q2)) + SH_MAG(1)*(P(21,3) + P(1,3)*SH_MAG(3) + P(2,3)*SH_MAG(2) + P(3,3)*SH_MAG(1) - P(18,3)*(SH_MAG(4) - SH_MAG(5) + SH_MAG(6) - SH_MAG(7)) - P(17,3)*(2*q0*q3 - 2*q1*q2) + P(19,3)*(2*q0*q1 + 2*q2*q3) - P(4,3)*(SH_MAG(8) + SH_MAG(9) - 2*magD*q2)) - (SH_MAG(4) - SH_MAG(5) + SH_MAG(6) - SH_MAG(7))*(P(21,18) + P(1,18)*SH_MAG(3) + P(2,18)*SH_MAG(2) + P(3,18)*SH_MAG(1) - P(18,18)*(SH_MAG(4) - SH_MAG(5) + SH_MAG(6) - SH_MAG(7)) - P(17,18)*(2*q0*q3 - 2*q1*q2) + P(19,18)*(2*q0*q1 + 2*q2*q3) - P(4,18)*(SH_MAG(8) + SH_MAG(9) - 2*magD*q2)) - P(4,21)*(SH_MAG(8) + SH_MAG(9) - 2*magD*q2));
        SK_MY(2) = SH_MAG(4) - SH_MAG(5) + SH_MAG(6) - SH_MAG(7);
        SK_MY(3) = SH_MAG(8) + SH_MAG(9) - 2*magD*q2;
        SK_MY(4) = 2*q0*q3 - 2*q1*q2;
        SK_MY(5) = 2*q0*q1 + 2*q2*q3;
        
        
    
        Kfusion(1) = SK_MY(1)*(P(1,21) + P(1,1)*SH_MAG(3) + P(1,2)*SH_MAG(2) + P(1,3)*SH_MAG(1) - P(1,4)*SK_MY(3) - P(1,18)*SK_MY(2) - P(1,17)*SK_MY(4) + P(1,19)*SK_MY(5));
        Kfusion(2) = SK_MY(1)*(P(2,21) + P(2,1)*SH_MAG(3) + P(2,2)*SH_MAG(2) + P(2,3)*SH_MAG(1) - P(2,4)*SK_MY(3) - P(2,18)*SK_MY(2) - P(2,17)*SK_MY(4) + P(2,19)*SK_MY(5));
        Kfusion(3) = SK_MY(1)*(P(3,21) + P(3,1)*SH_MAG(3) + P(3,2)*SH_MAG(2) + P(3,3)*SH_MAG(1) - P(3,4)*SK_MY(3) - P(3,18)*SK_MY(2) - P(3,17)*SK_MY(4) + P(3,19)*SK_MY(5));
        Kfusion(4) = SK_MY(1)*(P(4,21) + P(4,1)*SH_MAG(3) + P(4,2)*SH_MAG(2) + P(4,3)*SH_MAG(1) - P(4,4)*SK_MY(3) - P(4,18)*SK_MY(2) - P(4,17)*SK_MY(4) + P(4,19)*SK_MY(5));
        Kfusion(5) = SK_MY(1)*(P(5,21) + P(5,1)*SH_MAG(3) + P(5,2)*SH_MAG(2) + P(5,3)*SH_MAG(1) - P(5,4)*SK_MY(3) - P(5,18)*SK_MY(2) - P(5,17)*SK_MY(4) + P(5,19)*SK_MY(5));
        Kfusion(6) = SK_MY(1)*(P(6,21) + P(6,1)*SH_MAG(3) + P(6,2)*SH_MAG(2) + P(6,3)*SH_MAG(1) - P(6,4)*SK_MY(3) - P(6,18)*SK_MY(2) - P(6,17)*SK_MY(4) + P(6,19)*SK_MY(5));
        Kfusion(7) = SK_MY(1)*(P(7,21) + P(7,1)*SH_MAG(3) + P(7,2)*SH_MAG(2) + P(7,3)*SH_MAG(1) - P(7,4)*SK_MY(3) - P(7,18)*SK_MY(2) - P(7,17)*SK_MY(4) + P(7,19)*SK_MY(5));
        Kfusion(8) = SK_MY(1)*(P(8,21) + P(8,1)*SH_MAG(3) + P(8,2)*SH_MAG(2) + P(8,3)*SH_MAG(1) - P(8,4)*SK_MY(3) - P(8,18)*SK_MY(2) - P(8,17)*SK_MY(4) + P(8,19)*SK_MY(5));
        Kfusion(9) = SK_MY(1)*(P(9,21) + P(9,1)*SH_MAG(3) + P(9,2)*SH_MAG(2) + P(9,3)*SH_MAG(1) - P(9,4)*SK_MY(3) - P(9,18)*SK_MY(2) - P(9,17)*SK_MY(4) + P(9,19)*SK_MY(5));
        Kfusion(10) = SK_MY(1)*(P(10,21) + P(10,1)*SH_MAG(3) + P(10,2)*SH_MAG(2) + P(10,3)*SH_MAG(1) - P(10,4)*SK_MY(3) - P(10,18)*SK_MY(2) - P(10,17)*SK_MY(4) + P(10,19)*SK_MY(5));
        Kfusion(11) = SK_MY(1)*(P(11,21) + P(11,1)*SH_MAG(3) + P(11,2)*SH_MAG(2) + P(11,3)*SH_MAG(1) - P(11,4)*SK_MY(3) - P(11,18)*SK_MY(2) - P(11,17)*SK_MY(4) + P(11,19)*SK_MY(5));
        Kfusion(12) = SK_MY(1)*(P(12,21) + P(12,1)*SH_MAG(3) + P(12,2)*SH_MAG(2) + P(12,3)*SH_MAG(1) - P(12,4)*SK_MY(3) - P(12,18)*SK_MY(2) - P(12,17)*SK_MY(4) + P(12,19)*SK_MY(5));
        Kfusion(13) = SK_MY(1)*(P(13,21) + P(13,1)*SH_MAG(3) + P(13,2)*SH_MAG(2) + P(13,3)*SH_MAG(1) - P(13,4)*SK_MY(3) - P(13,18)*SK_MY(2) - P(13,17)*SK_MY(4) + P(13,19)*SK_MY(5));
        Kfusion(14) = SK_MY(1)*(P(14,21) + P(14,1)*SH_MAG(3) + P(14,2)*SH_MAG(2) + P(14,3)*SH_MAG(1) - P(14,4)*SK_MY(3) - P(14,18)*SK_MY(2) - P(14,17)*SK_MY(4) + P(14,19)*SK_MY(5));
        Kfusion(15) = SK_MY(1)*(P(15,21) + P(15,1)*SH_MAG(3) + P(15,2)*SH_MAG(2) + P(15,3)*SH_MAG(1) - P(15,4)*SK_MY(3) - P(15,18)*SK_MY(2) - P(15,17)*SK_MY(4) + P(15,19)*SK_MY(5));
        Kfusion(16) = SK_MY(1)*(P(16,21) + P(16,1)*SH_MAG(3) + P(16,2)*SH_MAG(2) + P(16,3)*SH_MAG(1) - P(16,4)*SK_MY(3) - P(16,18)*SK_MY(2) - P(16,17)*SK_MY(4) + P(16,19)*SK_MY(5));
        Kfusion(17) = SK_MY(1)*(P(17,21) + P(17,1)*SH_MAG(3) + P(17,2)*SH_MAG(2) + P(17,3)*SH_MAG(1) - P(17,4)*SK_MY(3) - P(17,18)*SK_MY(2) - P(17,17)*SK_MY(4) + P(17,19)*SK_MY(5));
        Kfusion(18) = SK_MY(1)*(P(18,21) + P(18,1)*SH_MAG(3) + P(18,2)*SH_MAG(2) + P(18,3)*SH_MAG(1) - P(18,4)*SK_MY(3) - P(18,18)*SK_MY(2) - P(18,17)*SK_MY(4) + P(18,19)*SK_MY(5));
        Kfusion(19) = SK_MY(1)*(P(19,21) + P(19,1)*SH_MAG(3) + P(19,2)*SH_MAG(2) + P(19,3)*SH_MAG(1) - P(19,4)*SK_MY(3) - P(19,18)*SK_MY(2) - P(19,17)*SK_MY(4) + P(19,19)*SK_MY(5));
        Kfusion(20) = SK_MY(1)*(P(20,21) + P(20,1)*SH_MAG(3) + P(20,2)*SH_MAG(2) + P(20,3)*SH_MAG(1) - P(20,4)*SK_MY(3) - P(20,18)*SK_MY(2) - P(20,17)*SK_MY(4) + P(20,19)*SK_MY(5));
        Kfusion(21) = SK_MY(1)*(P(21,21) + P(21,1)*SH_MAG(3) + P(21,2)*SH_MAG(2) + P(21,3)*SH_MAG(1) - P(21,4)*SK_MY(3) - P(21,18)*SK_MY(2) - P(21,17)*SK_MY(4) + P(21,19)*SK_MY(5));
        Kfusion(22) = SK_MY(1)*(P(22,21) + P(22,1)*SH_MAG(3) + P(22,2)*SH_MAG(2) + P(22,3)*SH_MAG(1) - P(22,4)*SK_MY(3) - P(22,18)*SK_MY(2) - P(22,17)*SK_MY(4) + P(22,19)*SK_MY(5));
        varInnov(2) = 1/SK_MY(1);
    elseif obsIndex == 3 % we ar now fusing the Z measurement
        % calculate the observation jacobian
        H_MAG(1) = SH_MAG(2);
        H_MAG(2) = -SH_MAG(3);
        H_MAG(3) = SH_MAG(8) + SH_MAG(9) - 2*magD*q2;
        H_MAG(4) = SH_MAG(1);
        H_MAG(17) = 2*q0*q2 + 2*q1*q3;
        H_MAG(18) = 2*q2*q3 - 2*q0*q1;
        H_MAG(19) = SH_MAG(4) - SH_MAG(5) - SH_MAG(6) + SH_MAG(7);
        H_MAG(22) = 1;
        % calculate the Kalman gain
        SK_MZ = zeros(5,1);
        SK_MZ(1) = 1/(P(22,22) + R_MAG + P(1,22)*SH_MAG(2) - P(2,22)*SH_MAG(3) + P(4,22)*SH_MAG(1) + P(19,22)*(SH_MAG(4) - SH_MAG(5) - SH_MAG(6) + SH_MAG(7)) + (2*q0*q2 + 2*q1*q3)*(P(22,17) + P(1,17)*SH_MAG(2) - P(2,17)*SH_MAG(3) + P(4,17)*SH_MAG(1) + P(19,17)*(SH_MAG(4) - SH_MAG(5) - SH_MAG(6) + SH_MAG(7)) + P(17,17)*(2*q0*q2 + 2*q1*q3) - P(18,17)*(2*q0*q1 - 2*q2*q3) + P(3,17)*(SH_MAG(8) + SH_MAG(9) - 2*magD*q2)) - (2*q0*q1 - 2*q2*q3)*(P(22,18) + P(1,18)*SH_MAG(2) - P(2,18)*SH_MAG(3) + P(4,18)*SH_MAG(1) + P(19,18)*(SH_MAG(4) - SH_MAG(5) - SH_MAG(6) + SH_MAG(7)) + P(17,18)*(2*q0*q2 + 2*q1*q3) - P(18,18)*(2*q0*q1 - 2*q2*q3) + P(3,18)*(SH_MAG(8) + SH_MAG(9) - 2*magD*q2)) + (SH_MAG(8) + SH_MAG(9) - 2*magD*q2)*(P(22,3) + P(1,3)*SH_MAG(2) - P(2,3)*SH_MAG(3) + P(4,3)*SH_MAG(1) + P(19,3)*(SH_MAG(4) - SH_MAG(5) - SH_MAG(6) + SH_MAG(7)) + P(17,3)*(2*q0*q2 + 2*q1*q3) - P(18,3)*(2*q0*q1 - 2*q2*q3) + P(3,3)*(SH_MAG(8) + SH_MAG(9) - 2*magD*q2)) + P(17,22)*(2*q0*q2 + 2*q1*q3) - P(18,22)*(2*q0*q1 - 2*q2*q3) + SH_MAG(2)*(P(22,1) + P(1,1)*SH_MAG(2) - P(2,1)*SH_MAG(3) + P(4,1)*SH_MAG(1) + P(19,1)*(SH_MAG(4) - SH_MAG(5) - SH_MAG(6) + SH_MAG(7)) + P(17,1)*(2*q0*q2 + 2*q1*q3) - P(18,1)*(2*q0*q1 - 2*q2*q3) + P(3,1)*(SH_MAG(8) + SH_MAG(9) - 2*magD*q2)) - SH_MAG(3)*(P(22,2) + P(1,2)*SH_MAG(2) - P(2,2)*SH_MAG(3) + P(4,2)*SH_MAG(1) + P(19,2)*(SH_MAG(4) - SH_MAG(5) - SH_MAG(6) + SH_MAG(7)) + P(17,2)*(2*q0*q2 + 2*q1*q3) - P(18,2)*(2*q0*q1 - 2*q2*q3) + P(3,2)*(SH_MAG(8) + SH_MAG(9) - 2*magD*q2)) + SH_MAG(1)*(P(22,4) + P(1,4)*SH_MAG(2) - P(2,4)*SH_MAG(3) + P(4,4)*SH_MAG(1) + P(19,4)*(SH_MAG(4) - SH_MAG(5) - SH_MAG(6) + SH_MAG(7)) + P(17,4)*(2*q0*q2 + 2*q1*q3) - P(18,4)*(2*q0*q1 - 2*q2*q3) + P(3,4)*(SH_MAG(8) + SH_MAG(9) - 2*magD*q2)) + (SH_MAG(4) - SH_MAG(5) - SH_MAG(6) + SH_MAG(7))*(P(22,19) + P(1,19)*SH_MAG(2) - P(2,19)*SH_MAG(3) + P(4,19)*SH_MAG(1) + P(19,19)*(SH_MAG(4) - SH_MAG(5) - SH_MAG(6) + SH_MAG(7)) + P(17,19)*(2*q0*q2 + 2*q1*q3) - P(18,19)*(2*q0*q1 - 2*q2*q3) + P(3,19)*(SH_MAG(8) + SH_MAG(9) - 2*magD*q2)) + P(3,22)*(SH_MAG(8) + SH_MAG(9) - 2*magD*q2));
        SK_MZ(2) = SH_MAG(4) - SH_MAG(5) - SH_MAG(6) + SH_MAG(7);
        SK_MZ(3) = SH_MAG(8) + SH_MAG(9) - 2*magD*q2;
        SK_MZ(4) = 2*q0*q1 - 2*q2*q3;
        SK_MZ(5) = 2*q0*q2 + 2*q1*q3;
        
   
        Kfusion(1) = SK_MZ(1)*(P(1,22) + P(1,1)*SH_MAG(2) - P(1,2)*SH_MAG(3) + P(1,4)*SH_MAG(1) + P(1,3)*SK_MZ(3) + P(1,19)*SK_MZ(2) + P(1,17)*SK_MZ(5) - P(1,18)*SK_MZ(4));
        Kfusion(2) = SK_MZ(1)*(P(2,22) + P(2,1)*SH_MAG(2) - P(2,2)*SH_MAG(3) + P(2,4)*SH_MAG(1) + P(2,3)*SK_MZ(3) + P(2,19)*SK_MZ(2) + P(2,17)*SK_MZ(5) - P(2,18)*SK_MZ(4));
        Kfusion(3) = SK_MZ(1)*(P(3,22) + P(3,1)*SH_MAG(2) - P(3,2)*SH_MAG(3) + P(3,4)*SH_MAG(1) + P(3,3)*SK_MZ(3) + P(3,19)*SK_MZ(2) + P(3,17)*SK_MZ(5) - P(3,18)*SK_MZ(4));
        Kfusion(4) = SK_MZ(1)*(P(4,22) + P(4,1)*SH_MAG(2) - P(4,2)*SH_MAG(3) + P(4,4)*SH_MAG(1) + P(4,3)*SK_MZ(3) + P(4,19)*SK_MZ(2) + P(4,17)*SK_MZ(5) - P(4,18)*SK_MZ(4));
        Kfusion(5) = SK_MZ(1)*(P(5,22) + P(5,1)*SH_MAG(2) - P(5,2)*SH_MAG(3) + P(5,4)*SH_MAG(1) + P(5,3)*SK_MZ(3) + P(5,19)*SK_MZ(2) + P(5,17)*SK_MZ(5) - P(5,18)*SK_MZ(4));
        Kfusion(6) = SK_MZ(1)*(P(6,22) + P(6,1)*SH_MAG(2) - P(6,2)*SH_MAG(3) + P(6,4)*SH_MAG(1) + P(6,3)*SK_MZ(3) + P(6,19)*SK_MZ(2) + P(6,17)*SK_MZ(5) - P(6,18)*SK_MZ(4));
        Kfusion(7) = SK_MZ(1)*(P(7,22) + P(7,1)*SH_MAG(2) - P(7,2)*SH_MAG(3) + P(7,4)*SH_MAG(1) + P(7,3)*SK_MZ(3) + P(7,19)*SK_MZ(2) + P(7,17)*SK_MZ(5) - P(7,18)*SK_MZ(4));
        Kfusion(8) = SK_MZ(1)*(P(8,22) + P(8,1)*SH_MAG(2) - P(8,2)*SH_MAG(3) + P(8,4)*SH_MAG(1) + P(8,3)*SK_MZ(3) + P(8,19)*SK_MZ(2) + P(8,17)*SK_MZ(5) - P(8,18)*SK_MZ(4));
        Kfusion(9) = SK_MZ(1)*(P(9,22) + P(9,1)*SH_MAG(2) - P(9,2)*SH_MAG(3) + P(9,4)*SH_MAG(1) + P(9,3)*SK_MZ(3) + P(9,19)*SK_MZ(2) + P(9,17)*SK_MZ(5) - P(9,18)*SK_MZ(4));
        Kfusion(10) = SK_MZ(1)*(P(10,22) + P(10,1)*SH_MAG(2) - P(10,2)*SH_MAG(3) + P(10,4)*SH_MAG(1) + P(10,3)*SK_MZ(3) + P(10,19)*SK_MZ(2) + P(10,17)*SK_MZ(5) - P(10,18)*SK_MZ(4));
        Kfusion(11) = SK_MZ(1)*(P(11,22) + P(11,1)*SH_MAG(2) - P(11,2)*SH_MAG(3) + P(11,4)*SH_MAG(1) + P(11,3)*SK_MZ(3) + P(11,19)*SK_MZ(2) + P(11,17)*SK_MZ(5) - P(11,18)*SK_MZ(4));
        Kfusion(12) = SK_MZ(1)*(P(12,22) + P(12,1)*SH_MAG(2) - P(12,2)*SH_MAG(3) + P(12,4)*SH_MAG(1) + P(12,3)*SK_MZ(3) + P(12,19)*SK_MZ(2) + P(12,17)*SK_MZ(5) - P(12,18)*SK_MZ(4));
        Kfusion(13) = SK_MZ(1)*(P(13,22) + P(13,1)*SH_MAG(2) - P(13,2)*SH_MAG(3) + P(13,4)*SH_MAG(1) + P(13,3)*SK_MZ(3) + P(13,19)*SK_MZ(2) + P(13,17)*SK_MZ(5) - P(13,18)*SK_MZ(4));
        Kfusion(14) = SK_MZ(1)*(P(14,22) + P(14,1)*SH_MAG(2) - P(14,2)*SH_MAG(3) + P(14,4)*SH_MAG(1) + P(14,3)*SK_MZ(3) + P(14,19)*SK_MZ(2) + P(14,17)*SK_MZ(5) - P(14,18)*SK_MZ(4));
        Kfusion(15) = SK_MZ(1)*(P(15,22) + P(15,1)*SH_MAG(2) - P(15,2)*SH_MAG(3) + P(15,4)*SH_MAG(1) + P(15,3)*SK_MZ(3) + P(15,19)*SK_MZ(2) + P(15,17)*SK_MZ(5) - P(15,18)*SK_MZ(4));
        Kfusion(16) = SK_MZ(1)*(P(16,22) + P(16,1)*SH_MAG(2) - P(16,2)*SH_MAG(3) + P(16,4)*SH_MAG(1) + P(16,3)*SK_MZ(3) + P(16,19)*SK_MZ(2) + P(16,17)*SK_MZ(5) - P(16,18)*SK_MZ(4));
        Kfusion(17) = SK_MZ(1)*(P(17,22) + P(17,1)*SH_MAG(2) - P(17,2)*SH_MAG(3) + P(17,4)*SH_MAG(1) + P(17,3)*SK_MZ(3) + P(17,19)*SK_MZ(2) + P(17,17)*SK_MZ(5) - P(17,18)*SK_MZ(4));
        Kfusion(18) = SK_MZ(1)*(P(18,22) + P(18,1)*SH_MAG(2) - P(18,2)*SH_MAG(3) + P(18,4)*SH_MAG(1) + P(18,3)*SK_MZ(3) + P(18,19)*SK_MZ(2) + P(18,17)*SK_MZ(5) - P(18,18)*SK_MZ(4));
        Kfusion(19) = SK_MZ(1)*(P(19,22) + P(19,1)*SH_MAG(2) - P(19,2)*SH_MAG(3) + P(19,4)*SH_MAG(1) + P(19,3)*SK_MZ(3) + P(19,19)*SK_MZ(2) + P(19,17)*SK_MZ(5) - P(19,18)*SK_MZ(4));
        Kfusion(20) = SK_MZ(1)*(P(20,22) + P(20,1)*SH_MAG(2) - P(20,2)*SH_MAG(3) + P(20,4)*SH_MAG(1) + P(20,3)*SK_MZ(3) + P(20,19)*SK_MZ(2) + P(20,17)*SK_MZ(5) - P(20,18)*SK_MZ(4));
        Kfusion(21) = SK_MZ(1)*(P(21,22) + P(21,1)*SH_MAG(2) - P(21,2)*SH_MAG(3) + P(21,4)*SH_MAG(1) + P(21,3)*SK_MZ(3) + P(21,19)*SK_MZ(2) + P(21,17)*SK_MZ(5) - P(21,18)*SK_MZ(4));
        Kfusion(22) = SK_MZ(1)*(P(22,22) + P(22,1)*SH_MAG(2) - P(22,2)*SH_MAG(3) + P(22,4)*SH_MAG(1) + P(22,3)*SK_MZ(3) + P(22,19)*SK_MZ(2) + P(22,17)*SK_MZ(5) - P(22,18)*SK_MZ(4));
        varInnov(3) = 1/SK_MZ(1);
    end
    % Calculate the measurement innovation
    innovation(obsIndex) = MagPred(obsIndex) - MagData(obsIndex);
    % Check the innovation for consistencey and don't fuse if > 5Sigma
    if ((innovation(obsIndex)^2) / varInnov(obsIndex)) < 25.0
        xk = Kfusion * innovation(obsIndex);
        states = states - xk;
        % normalise the quaternion states
        quatMag = sqrt(states(1)^2 + states(2)^2 + states(3)^2 + states(4)^2);
        if (quatMag > 1e-12)
            states(1:4) = states(1:4) / quatMag;
        end
        % correct the covariance P = (I - K*H)*P
        % take advantage of the empty columns in KH to reduce the
        % number of operations
        for i = 1:22
            for j = 1:4
                KH(i,j) = Kfusion(i,1)*H_MAG(1,j);
            end
            for j = 17:22
                KH(i,j) = Kfusion(i,1)*H_MAG(1,j);
            end
        end
        for i = 1:22
            for j = 1:22
                for k = 1:4
                    KHP(i,j) = KHP(i,j) + KH(i,k)*P(k,j);
                end
                for k = 17:22
                    KHP(i,j) = KHP(i,j) + KH(i,k)*P(k,j);
                end
            end
        end
        P = P - KHP;
    end
    % increment observation index
    obsIndex = obsIndex + 1;
end

end