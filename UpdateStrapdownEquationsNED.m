%#codegen
function [states,correctedDelAng,correctedDelVel,accNavMag]  = UpdateStrapdownEquationsNED( ...
    states, ...
    dAngIMU, ...
    dVelIMU, ...
    dtIMU)
% define persistent variables for previous delta angle and velocity which
% are required for sculling and coning error corrections
persistent prevDelAng;
if isempty(prevDelAng)
    prevDelAng = single([0,0,0]);
end
persistent prevDelVel;
if isempty(prevDelVel)
    prevDelVel = single([0,0,0]);
end
persistent Tbn;
if isempty(Tbn)
    Tbn = single(eye(3));
end

gravityNED = [0.0;0.0;9.80665];

% Remove sensor bias errors
correctedDelAng = double(dAngIMU) - states(11:13)';
correctedDelVel = double(dVelIMU) - states(14:16)';


% Apply corrections for earths rotation rate and coning errors
 correctedDelAng   = correctedDelAng  + 8.333333333333333e-2*cross(prevDelAng , correctedDelAng );
% Apply rotational and skulling corrections
correctedDelVel= correctedDelVel+ ...
            0.5*cross( correctedDelAng , correctedDelVel)+...% rotational correction
            2/3*(cross(prevDelVel  , correctedDelAng) - cross(correctedDelVel ,prevDelAng)); % sculling correction
        
% %  Apply rotational and skulling corrections
%          correctedDelVel= dVel+ ...
%             0.5*cross(prevDelAng + dAng , prevDelVel + dVel)+ 1/6*cross(prevDelAng + correctedDelAng , cross(prevDelAng + correctedDelAng , prevDelVel + correctedDelVel)) + ... % rotational correction
%             1/12*(cross(prevDelAng , correctedDelVel) + cross(prevDelVel , correctedDelAng)); % sculling correction,1*3


% Save current measurements
prevDelAng = correctedDelAng;
prevDelVel = correctedDelVel;
% Convert the rotation vector to its equivalent quaternion
rotationMag = sqrt(correctedDelAng(1)^2 + correctedDelAng(2)^2 + correctedDelAng(3)^2);
if rotationMag<1e-12
  deltaQuat = single([1,0,0,0]);
else
  deltaQuat = [cos(0.5*rotationMag), correctedDelAng/rotationMag*sin(0.5*rotationMag)]';
end

% Update the quaternions by rotating from the previous attitude through
% the delta angle rotation quaternion
 states(1:4,1)= [states(1)*deltaQuat(1)-transpose(states(2:4))*deltaQuat(2:4); states(1)*deltaQuat(2:4) + deltaQuat(1)*states(2:4) + cross(states(2:4),deltaQuat(2:4))];

% Normalise the quaternions and update the quaternion states
quatMag = sqrt(states(1,1)^2 + states(2,1)^2 + states(3,1)^2 + states(4,1)^2);
if (quatMag < 1e-16)
    states(1:4) =  states(1:4);
else
    states(1:4) =  states(1:4) / quatMag;
end

% Calculate the body to nav cosine matrix
Tbn = [states(1)^2 + states(2)^2 - states(3)^2 - states(4)^2, 2*(states(2)*states(3) - states(1)*states(4)), 2*(states(2)*states(4) + states(1)*states(3));...
      2*(states(2)*states(3) + states(1)*states(4)), states(1)^2 - states(2)^2 + states(3)^2 - states(4)^2, 2*(states(3)*states(4) - states(1)*states(2));...
      2*(states(2)*states(4)-states(1)*states(3)), 2*(states(3)*states(4) + states(1)*states(2)), states(1)^2 - states(2)^2 - states(3)^2 + states(4)^2];
  
% transform body delta velocities to delta velocities in the nav frame
correctedDelVel =correctedDelVel';
delVelNav = Tbn * correctedDelVel + gravityNED*dtIMU;

% calculate the magnitude of the nav acceleration (required for GPS
% variance estimation)
accNavMag = sqrt(delVelNav(1)^2 + delVelNav(2)^2 + delVelNav(3)^2) / dtIMU;

% If calculating position save previous velocity
lastVelocity = states(5:7);

% Sum delta velocities to get velocity after removing gravitational skew and correcting for transport rate
states(5:7) = states(5:7) + delVelNav;


% If calculating postions, do a trapezoidal integration for position
states(8:10)     = states(8:10) + 0.5*(states(5:7) + lastVelocity)*dtIMU;

% Copy the remaining states across
states(11:16) = states(11:16);

end