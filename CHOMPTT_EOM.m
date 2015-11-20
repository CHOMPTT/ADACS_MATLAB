% CHOMPTT_EOM.m -----------------------------------------------------------
% 
% This file is used to integrate the satellite dynamics
%
% INPUT:
%   t is the current time
%   t_step is the time step
%   state is the state vector (n x 1)
%   Tc is the control torque (3 x 1)
%   Fd is the disturbance vector (3 x 1)
%   Td is the disturbance torque (3 x 1)
%
% OUTPUT:
%   state derivative
%

function [state_dot] = CHOMPTT_EOM( t, state, Tc, Fd, Td, w_i_e);     

% Get State Variables
x  =  state(1);        
y  =  state(2);        
z  =  state(3);         
vel_x =  state(4);         
vel_y =  state(5);         
vel_z =  state(6);        
omega_x =  state(7);        
omega_y =  state(8);         
omega_z =  state(9);         
q1 = state(10);        
q2 = state(11);        
q3 = state(12);         
q4 = state(13);         
r_sc = state(1:3);
wx = omega_x;
wy=omega_y;
wz=omega_z;


q = [q1;q2;q3;q4];
omega = [omega_x; omega_y;omega_z];
% Constants
% List any constants you might use here 
R_s_e = [q(1)^2-q(2)^2-q(3)^2+q(4)^2, 2*(q(1)*q(2)+q(3)*q(4)), 2*(q(1)*q(3)-q(2)*q(4)); 2*(q(1)*q(2)-q(3)*q(4)), -q(1)^2+q(2)^2-q(3)^2+q(4)^2, 2*(q(2)*q(3)+q(1)*q(4)); 2*(q(1)*q(3)+q(2)*q(4)), 2*(q(2)*q(3)+q(1)*q(4)), -q(1)^2-q(2)^2+q(3)^2+q(4)^4] ;
GM = 398600.4415 *(10^9); %m^3/s^2


% Sat inertia constants
Msc=3;     %Spacecraft Mass := 3Kg
a=0.1;      %Spacecraft length along x-axis  (m)
b=0.1;      %Spacecraft length along y-axis  (m)
c=0.34;    %Spacecraft length along z-axis  (m)

%Spacecraft Moments of Inertia
Ix=1/12*Msc*(b^2+c^2) ;            
Iy=1/12*Msc*(a^2+c^2) ;
Iz=1/12*Msc*(a^2+b^2) ;
inertia_matrix = [Ix 0 0; 0 Iy 0; 0 0 Iz];

% Control forces
Tcx = Tc(1);
Tcy = Tc(2);
Tcz = Tc(3);

% Disturbance forces and torques
Fdx = Fd(1);    Tdx = Td(1);
Fdy = Fd(2);    Tdy = Td(2);
Fdz = Fd(3);    Tdz = Td(3);


% Quaternion Derivatives (EOM given)
q_dot = 0.5*cross(omega,[q1;q2;q3]) + 0.5*q(4)*omega;
q1_dot = q_dot(1) ;
q2_dot = q_dot(2);
q3_dot = q_dot(3);
q4_dot = -0.5*dot(omega,[q1;q2;q3]);

% Angular Acceleration from Rotational EOM
omega_dot = inv(inertia_matrix)*(cross((-1*(R_s_e*w_i_e+omega)),inertia_matrix*(R_s_e*w_i_e+omega))-(Tc+Td));
omega_dot_x = -(Tcx + Tdx - Iy*wy(t)*wz(t) + Iz*wy(t)*wz(t) + Iy*We*q1(t)^2*wy(t) + Iy*We*q2(t)^2*wy(t) - Iy*We*q3(t)^2*wy(t) - Iy*We*q4(t)^4*wy(t) - Iz*We*q1(t)^2*wy(t) - Iz*We*q2(t)^2*wy(t) + Iz*We*q3(t)^2*wy(t) + Iz*We*q4(t)^4*wy(t) + 2*Ix*We*dq1*q3(t) + 2*Ix*We*dq3*q1(t) - 2*Ix*We*dq2*q4(t) - 2*Ix*We*dq4*q2(t) - 2*Iy*We^2*q2(t)*q3(t)^3 + 2*Iy*We^2*q1(t)^3*q4(t) + 2*Iy*We^2*q2(t)^3*q3(t) - 2*Iy*We^2*q1(t)*q4(t)^5 + 2*Iz*We^2*q2(t)*q3(t)^3 - 2*Iz*We^2*q1(t)^3*q4(t) - 2*Iz*We^2*q2(t)^3*q3(t) + 2*Iz*We^2*q1(t)*q4(t)^5 + 2*Iy*We^2*q1(t)^2*q2(t)*q3(t) + 2*Iy*We^2*q1(t)*q2(t)^2*q4(t) - 2*Iy*We^2*q1(t)*q3(t)^2*q4(t) - 2*Iy*We^2*q2(t)*q3(t)*q4(t)^4 - 2*Iz*We^2*q1(t)^2*q2(t)*q3(t) - 2*Iz*We^2*q1(t)*q2(t)^2*q4(t) + 2*Iz*We^2*q1(t)*q3(t)^2*q4(t) + 2*Iz*We^2*q2(t)*q3(t)*q4(t)^4 - 2*Iy*We*q1(t)*q4(t)*wz(t) - 2*Iy*We*q2(t)*q3(t)*wz(t) + 2*Iz*We*q1(t)*q4(t)*wz(t) + 2*Iz*We*q2(t)*q3(t)*wz(t))/Ix;
omega_dot_y = -(Tcy + Tdy + Ix*wx(t)*wz(t) - Iz*wx(t)*wz(t) - Ix*We*q1(t)^2*wx(t) - Ix*We*q2(t)^2*wx(t) + Ix*We*q3(t)^2*wx(t) + Ix*We*q4(t)^4*wx(t) + Iz*We*q1(t)^2*wx(t) + Iz*We*q2(t)^2*wx(t) - Iz*We*q3(t)^2*wx(t) - Iz*We*q4(t)^4*wx(t) + 2*Iy*We*dq1*q4(t) + 2*Iy*We*dq2*q3(t) + 2*Iy*We*dq3*q2(t) + 2*Iy*We*dq4*q1(t) + 2*Ix*We^2*q1(t)*q3(t)^3 - 2*Ix*We^2*q1(t)^3*q3(t) + 2*Ix*We^2*q2(t)^3*q4(t) - 2*Ix*We^2*q2(t)*q4(t)^5 - 2*Iz*We^2*q1(t)*q3(t)^3 + 2*Iz*We^2*q1(t)^3*q3(t) - 2*Iz*We^2*q2(t)^3*q4(t) + 2*Iz*We^2*q2(t)*q4(t)^5 - 2*Ix*We^2*q1(t)*q2(t)^2*q3(t) + 2*Ix*We^2*q1(t)^2*q2(t)*q4(t) - 2*Ix*We^2*q2(t)*q3(t)^2*q4(t) + 2*Ix*We^2*q1(t)*q3(t)*q4(t)^4 + 2*Iz*We^2*q1(t)*q2(t)^2*q3(t) - 2*Iz*We^2*q1(t)^2*q2(t)*q4(t) + 2*Iz*We^2*q2(t)*q3(t)^2*q4(t) - 2*Iz*We^2*q1(t)*q3(t)*q4(t)^4 + 2*Ix*We*q1(t)*q3(t)*wz(t) - 2*Ix*We*q2(t)*q4(t)*wz(t) - 2*Iz*We*q1(t)*q3(t)*wz(t) + 2*Iz*We*q2(t)*q4(t)*wz(t))/Iy;
omega_dot_z = -(Tcz + Tdz - Ix*wx(t)*wy(t) + Iy*wx(t)*wy(t) - 2*Iz*We*dq1*q1(t) - 2*Iz*We*dq2*q2(t) + 2*Iz*We*dq3*q3(t) + 4*Iz*We*dq4*q4(t)^3 - 4*Ix*We^2*q1(t)*q2(t)*q3(t)^2 + 4*Ix*We^2*q1(t)*q2(t)*q4(t)^2 - 4*Ix*We^2*q1(t)^2*q3(t)*q4(t) + 4*Ix*We^2*q2(t)^2*q3(t)*q4(t) + 4*Iy*We^2*q1(t)*q2(t)*q3(t)^2 - 4*Iy*We^2*q1(t)*q2(t)*q4(t)^2 + 4*Iy*We^2*q1(t)^2*q3(t)*q4(t) - 4*Iy*We^2*q2(t)^2*q3(t)*q4(t) - 2*Ix*We*q1(t)*q4(t)*wx(t) - 2*Ix*We*q2(t)*q3(t)*wx(t) - 2*Ix*We*q1(t)*q3(t)*wy(t) + 2*Iy*We*q1(t)*q4(t)*wx(t) + 2*Iy*We*q2(t)*q3(t)*wx(t) + 2*Ix*We*q2(t)*q4(t)*wy(t) + 2*Iy*We*q1(t)*q3(t)*wy(t) - 2*Iy*We*q2(t)*q4(t)*wy(t))/Iz;


% Acceleration from Translational EOM
acc_x = Fdx/Msc - (GM/norm(r_sc)^3)*x;
acc_y = Fdy/Msc - (GM/norm(r_sc)^3)*y;
acc_z = Fdz/Msc - (GM/norm(r_sc)^3)*z;

% Build the vectors
vel   = [vel_x; vel_y; vel_z];
accel = [acc_x; acc_y; acc_z];
omega_dot = [omega_dot_x; omega_dot_y; omega_dot_z];
q_dot = [q1_dot; q2_dot; q3_dot; q4_dot];

% Return the derivative of the state 
state_dot = [vel; accel; omega_dot; q_dot];

end %CHOMPTT_EOM.m -------------------------------------------------------
