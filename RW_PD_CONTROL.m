% RW_PD_CONTROL.m -----------------------------------------------------------
% 
% This file is used to simulate the reaction wheel control
%
% INPUT:
%   k is the proportional gain (1 x 3) [Kp_phi, Kp_theta, Kp_psi]
%   c is the derivative gain (1 x 3) [Kd_phi, Kd_theta, Kd_psi]
%   state is the state vector (13 x 1)
    %   state(1)=r_sc(1)=x position of satellite in ECI
    %   state(2)=r_sc(2)=y position of satellite in ECI
    %   state(3)=r_sc(3)=z position of satellite in ECI
    %   state(4)=vel_sc(1)= x tangential velocity of satellite in ECI
    %   state(5)=vel_sc(2)= y tangential velocity of satellite in ECI
    %   state(6)=vel_sc(3)= z tangential velocity of satellite in ECI
    %   state(7)=w_e_s(1)= x S/C angular velocity w.r.t ECEF in S/C frame
    %   state(8)=w_e_s(2)= x S/C angular velocity w.r.t ECEF in S/C frame
    %   state(9)=w_e_s(3)= x S/C angular velocity w.r.t ECEF in S/C frame
    %   state(10)=q(1)
    %   state(11)=q(2)
    %   state(12)=q(3)
    %   state(13)=q(4)
%   MAX_TORQUE = scalar maximum torque of the wheel
%
% OUTPUT:
%   Tc is the control torque (3 x 1)
%

function Tc = RW_PD_CONTROL(k, c, state, MAX_TORQUE);
qt = transpose(state(10:13));
error = Quat2Eul321(qt); 
error_dot = state(7:9);
%max_torque = 3.75e-3;

T_x = -k(1)*(error(3)) - c(1)*(error_dot(1));
T_y = -k(2)*(error(2)) - c(2)*(error_dot(2));
T_z = -k(3)*(error(1)) - c(3)*(error_dot(3));
T = [T_x; T_y; T_z];

% Checking for Torque Max-Min Requirements

for j = 1:3
   if abs(T(j,1)) >= MAX_TORQUE
        if sign(T(j,1)) == 1
            T(j,1) = MAX_TORQUE;
        else
            T(j,1) = -MAX_TORQUE;
        end
    end
end
Tc = T;
end
