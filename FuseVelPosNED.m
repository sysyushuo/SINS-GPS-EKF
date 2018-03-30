%#codegen
function [...
    states, ... % state vector after fusion of measurements (q0, ... , q3, VelN,E,D, PosN,E,D, WindN,E, P0, camera roll pitch yaw misalignment)
    P, ... % state covariance matrix after fusion of corrections
    innovation,... % NED velocity and NED position innovations (m/s, m)
    varInnov] ... % NED velocity and NED position innovation variance ((m/s)^2, m^2)
    = FuseVelPosNED( ...
    states, ... % predicted states from the INS (q0, ... , q3, VelN,E,D, PosN,E,D, WindN,E, P0, camera roll pitch yaw misalignment)
    P, ... % predicted covariance
    accNavMag, ... % magnitude of rate of change of velocity in navigation frame (m/s^2)
    FuseVelData, ... % Boolean to initiate fusion of velocity measurements
    VelNED, ... % NED velocity measurements (m/s)
    FusePosData,  ... % Boolean to initiate fusion of position measurements
    PosNE, ... % NE position measurements (m)
    FuseHgtData, ... % Boolean to initiate fusion of height measurements
    HgtMea) % D position measurement (m)
  

% initialise measurement health flags to false
velHealth = 0;
posHealth = 0;
hgtHealth = 0;

% Specify GPS velocity induced error
velErr = single(0.15*accNavMag);
posErr = single(0.15*accNavMag);
% Specify the GPS and height measurement variances.
R_OBS = single([(0.1^2 + velErr^2) (0.1^2 + velErr^2) (0.1^2 + velErr^2) (2.0^2 + posErr^2) (2.0^2 + posErr^2) 2.0^2]);
% initialise the innovation vector
innovation = single(zeros(1,6));

% initialise the innovation variance vector
varInnov = single(zeros(1,6));

% Define the observation vector
observation = single([VelNED,PosNE,-HgtMea]);

% Perform sequential fusion of GPS and height measurements. This assumes that the
% errors in the dfferent velocity and position componenets are
% uncorrelated which is not true, however in the absence of covariance
% data from the GPS receiver it is the only assumption we can make
% so we might as well take advantage of the computational efficiencies
% associated with sequential fusion
if (FuseVelData || FusePosData || FuseHgtData)
    % calculate innovation variances
    for obsIndex = 1:6
        stateIndex = 4 + obsIndex;
        varInnov(obsIndex) = P(stateIndex,stateIndex) + R_OBS(obsIndex);
    end
    % calculate innovations and check validity against limits
    if FuseVelData
        velInnov = states(5:7) - VelNED';
        if (velInnov(1)*velInnov(1) + velInnov(2)*velInnov(2) + velInnov(3)*velInnov(3)) < 25.0*(varInnov(1) + varInnov(2) + varInnov(3))      
            velHealth = 1;   
        else
            velHealth = 0;
        end
    end
    if FusePosData
        posInnov = states(8:9) - PosNE';
        if (posInnov(1)*posInnov(1) + posInnov(2)*posInnov(2)) < 100.0*(varInnov(4) + varInnov(5))
            posHealth = 1; 
        else
            posHealth = 0;
        end
    end
    if  FuseHgtData
        hgtInnov = states(10) + HgtMea;
        if (hgtInnov*hgtInnov) < 25.0*varInnov(6) 
            hgtHealth = 1;
        else
            hgtHealth = 0;          
        end
    end
    % Set range for sequential fusion of velocity and position measurements
    fuseData = zeros(1,6);
    if FuseVelData
        fuseData(1) = 1;
        fuseData(2) = 1;
        fuseData(3) = 1;
    end
    if FusePosData
        fuseData(4) = 1;
        fuseData(5) = 1;
    end
    if FuseHgtData
        fuseData(6) = 1;
    end
    % Fuse measurements sequentially
    for obsIndex = 1:6
        if fuseData(obsIndex)
            % Apply data health checks
            % If no airspeed data, start using GPS data after x seconds
            if (velHealth && (obsIndex >= 1 && obsIndex <= 3)) || ...
                    (posHealth && (obsIndex == 4 || obsIndex == 5)) || ...
                    (hgtHealth && (obsIndex == 6))
                stateIndex = 4 + obsIndex;
                % Calculate the measurement innovation, using states from the
                % appropriate time coordinate
                if (obsIndex >= 1 && obsIndex <= 3)
                    innovation(obsIndex) = states(stateIndex) - observation(obsIndex);
                elseif (obsIndex == 4 || obsIndex == 5)
                    innovation(obsIndex) = states(stateIndex) - observation(obsIndex);
                elseif (obsIndex == 6)
                    innovation(obsIndex) = states(stateIndex) - observation(obsIndex);
                else
                    innovation(obsIndex) = 0.0;
                end
                
                % Calculate the Kalman Gain
                varInnov(obsIndex) = P(stateIndex,stateIndex) + R_OBS(obsIndex); % data logging only
                SK = 1/varInnov(obsIndex);
                K = P(:,stateIndex)*SK;
                % Calculate state corrections
                xk = K * innovation(obsIndex);
                % Update the covariance - take advantage of direct observation of a
                % single state at stateIndex
                %P = (I - K*H)*P;
                KHP = single(zeros(16,16));
                for rowIndex = 1:16
                    for colIndex = 1:16
                        KHP(rowIndex,colIndex) = K(rowIndex)*P(stateIndex,colIndex);
                    end
                end
                P = P - KHP;
                % Apply the state corrections and re-normaliase the quaternions
                states = states - xk;
                % Normalise the quaternion states
                quatMag = sqrt(states(1)^2 + states(2)^2 + states(3)^2 + states(4)^2);
                if (quatMag > 1e-12)
                    states(1:4) = states(1:4) / quatMag;
                end
            end
        end
    end
end
end