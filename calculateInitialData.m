function covariance = calculateInitialData()

% Define the initial Euler angle covariance (Phi, Theta, Psi)
 InitialEulerCovariance  = single([(1*pi/180); (1*pi/180); (1*pi/180)].^2);
% Define the transformation vector from a 321 sequence Euler rotation vector to a q0,...,q3 quaternion vector
% linearised around a level orientation (roll,pitch = 0,0)
J_eul2quat = ...
    single([[ 0.0, 0.0, 0.0]; ...
    [ 0.5, 0.0, 0.0]; ...
    [ 0.0, 0.5, 0.0]; ...
    [ 0.0, 0.0, 0.5]]);

% Form the covariance matrix for the intial Euler angle coordinates
angleCov = diag(InitialEulerCovariance);

% Transform the Euler angle covariances into the equivalent quaternion covariances
quatCov = J_eul2quat*angleCov*transpose(J_eul2quat);

% define the state covariances with the exception of the quaternion covariances
Sigma_vel_NE = single(0.7); % 1 sigma uncertainty in horizontal velocity components
Sigma_vel_D  = single(0.7); % 1 sigma uncertainty in vertical velocity
Sigma_pos_NE = single(15); % 1 sigma uncertainty in horizontal position components
Sigma_pos_D  = single(15); % 1 sigma uncertainty in vertical position
Sigma_dAng   = single(1/10*pi/180*0.02); % 1 Sigma uncertainty in delta angle bias
Sigma_dVel   = single(0.1*0.02); % 1 Sigma uncertainty in delta velocity bias
% Sigma_MagNED = single(20);
% Sigma_MagXYZ = single(20);
covariance   = single(diag([0;0;0;0;Sigma_vel_NE*[1;1];Sigma_vel_D;Sigma_pos_NE*[1;1];Sigma_pos_D;Sigma_dAng*[1;1;1];Sigma_dVel*[1;1;1]].^2));
% Add the quaternion covariances
 covariance(1:4,1:4) = quatCov;

end

