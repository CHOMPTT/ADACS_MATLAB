% CHOMPTT_EOM.m -----------------------------------------------------------
% This file is used to integrate the satellite dynamics

% INPUT:
%   t is the current time
%   t_step is the time step
%   state is the state vector (n x 1)
%   Tc is the control torque (3 x 1)
%   Fd is the disturbance vector (3 x 1)
%   Td is the disturbance torque (3 x 1)

% OUTPUT:
%   state derivative

function [state_dot] = CHOMPTT_EOM(state, Tc, Fd, Td, w_i_e)     

%GET STATE VARIABLES
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
wy = omega_y;
wz = omega_z;
q = [q1;q2;q3;q4];
omega = [omega_x; omega_y;omega_z];

%CONSTANT(s)
%List any constants you might use here 
R_s_e = [q(1)^2-q(2)^2-q(3)^2+q(4)^2, 2*(q(1)*q(2)+q(3)*q(4)), 2*(q(1)*q(3)-q(2)*q(4)); 2*(q(1)*q(2)-q(3)*q(4)), -q(1)^2+q(2)^2-q(3)^2+q(4)^2, 2*(q(2)*q(3)+q(1)*q(4)); 2*(q(1)*q(3)+q(2)*q(4)), 2*(q(2)*q(3)+q(1)*q(4)), -q(1)^2-q(2)^2+q(3)^2+q(4)^4] ;
GM = 398600.4415 *(10^9); %m^3/s^2

%SAT INERTIA CONSTANTS
Msc=3;%Spacecraft Mass := 3Kg
a=0.1;%Spacecraft length along x-axis  (m)
b=0.1;%Spacecraft length along y-axis  (m)
c=0.34;%Spacecraft length along z-axis  (m)

%SPACECRAFT MOMENTS OF INERTIA
Ix=1/12*Msc*(b^2+c^2) ;            
Iy=1/12*Msc*(a^2+c^2) ;
Iz=1/12*Msc*(a^2+b^2) ;
inertia_matrix = [Ix 0 0; 0 Iy 0; 0 0 Iz];

%CONTROL FORCES
Tcx = Tc(1);
Tcy = Tc(2);
Tcz = Tc(3);

%DISTURBANCE FORCES AND TORQUES
Fdx = Fd(1); 
Fdy = Fd(2);    
Fdz = Fd(3); 
Tdx = Td(1);
Tdy = Td(2);
Tdz = Td(3);

%QUATERNION DERIVATIVES (EOM GIVEN)
dq = 0.5*cross(omega,[q1;q2;q3]) + 0.5*q(4)*omega;
dq1 = dq(1) ;
dq2 = dq(2);
dq3 = dq(3);
dq4 = -0.5*dot(omega,[q1;q2;q3]);

%ANGULAR ACCELERATION FROM ROTATIONAL EOM
omega_dot = inv(inertia_matrix)*(cross((-1*(R_s_e*w_i_e+omega)),inertia_matrix*(R_s_e*w_i_e+omega))-(Tc+Td));
omega_dot_x = -(Tcx + Tdx - Iy*wy*wz + Iz*wy*wz + Iy*w_i_e(3)*q1^2*wy + Iy*w_i_e(3)*q2^2*wy - Iy*w_i_e(3)*q3^2*wy - Iy*w_i_e(3)*q4^4*wy - Iz*w_i_e(3)*q1^2*wy - Iz*w_i_e(3)*q2^2*wy + Iz*w_i_e(3)*q3^2*wy + Iz*w_i_e(3)*q4^4*wy + 2*Ix*w_i_e(3)*dq1*q3 + 2*Ix*w_i_e(3)*dq3*q1 - 2*Ix*w_i_e(3)*dq2*q4 - 2*Ix*w_i_e(3)*dq4*q2 - 2*Iy*w_i_e(3)^2*q2*q3^3 + 2*Iy*w_i_e(3)^2*q1^3*q4 + 2*Iy*w_i_e(3)^2*q2^3*q3 - 2*Iy*w_i_e(3)^2*q1*q4^5 + 2*Iz*w_i_e(3)^2*q2*q3^3 - 2*Iz*w_i_e(3)^2*q1^3*q4 - 2*Iz*w_i_e(3)^2*q2^3*q3 + 2*Iz*w_i_e(3)^2*q1*q4^5 + 2*Iy*w_i_e(3)^2*q1^2*q2*q3 + 2*Iy*w_i_e(3)^2*q1*q2^2*q4 - 2*Iy*w_i_e(3)^2*q1*q3^2*q4 - 2*Iy*w_i_e(3)^2*q2*q3*q4^4 - 2*Iz*w_i_e(3)^2*q1^2*q2*q3 - 2*Iz*w_i_e(3)^2*q1*q2^2*q4 + 2*Iz*w_i_e(3)^2*q1*q3^2*q4 + 2*Iz*w_i_e(3)^2*q2*q3*q4^4 - 2*Iy*w_i_e(3)*q1*q4*wz - 2*Iy*w_i_e(3)*q2*q3*wz + 2*Iz*w_i_e(3)*q1*q4*wz + 2*Iz*w_i_e(3)*q2*q3*wz)/Ix;
omega_dot_y = -(Tcy + Tdy + Ix*wx*wz - Iz*wx*wz - Ix*w_i_e(3)*q1^2*wx - Ix*w_i_e(3)*q2^2*wx + Ix*w_i_e(3)*q3^2*wx + Ix*w_i_e(3)*q4^4*wx + Iz*w_i_e(3)*q1^2*wx + Iz*w_i_e(3)*q2^2*wx - Iz*w_i_e(3)*q3^2*wx - Iz*w_i_e(3)*q4^4*wx + 2*Iy*w_i_e(3)*dq1*q4 + 2*Iy*w_i_e(3)*dq2*q3 + 2*Iy*w_i_e(3)*dq3*q2 + 2*Iy*w_i_e(3)*dq4*q1 + 2*Ix*w_i_e(3)^2*q1*q3^3 - 2*Ix*w_i_e(3)^2*q1^3*q3 + 2*Ix*w_i_e(3)^2*q2^3*q4 - 2*Ix*w_i_e(3)^2*q2*q4^5 - 2*Iz*w_i_e(3)^2*q1*q3^3 + 2*Iz*w_i_e(3)^2*q1^3*q3 - 2*Iz*w_i_e(3)^2*q2^3*q4 + 2*Iz*w_i_e(3)^2*q2*q4^5 - 2*Ix*w_i_e(3)^2*q1*q2^2*q3 + 2*Ix*w_i_e(3)^2*q1^2*q2*q4 - 2*Ix*w_i_e(3)^2*q2*q3^2*q4 + 2*Ix*w_i_e(3)^2*q1*q3*q4^4 + 2*Iz*w_i_e(3)^2*q1*q2^2*q3 - 2*Iz*w_i_e(3)^2*q1^2*q2*q4 + 2*Iz*w_i_e(3)^2*q2*q3^2*q4 - 2*Iz*w_i_e(3)^2*q1*q3*q4^4 + 2*Ix*w_i_e(3)*q1*q3*wz - 2*Ix*w_i_e(3)*q2*q4*wz - 2*Iz*w_i_e(3)*q1*q3*wz + 2*Iz*w_i_e(3)*q2*q4*wz)/Iy;
omega_dot_z = -(Tcz + Tdz - Ix*wx*wy + Iy*wx*wy - 2*Iz*w_i_e(3)*dq1*q1 - 2*Iz*w_i_e(3)*dq2*q2 + 2*Iz*w_i_e(3)*dq3*q3 + 4*Iz*w_i_e(3)*dq4*q4^3 - 4*Ix*w_i_e(3)^2*q1*q2*q3^2 + 4*Ix*w_i_e(3)^2*q1*q2*q4^2 - 4*Ix*w_i_e(3)^2*q1^2*q3*q4 + 4*Ix*w_i_e(3)^2*q2^2*q3*q4 + 4*Iy*w_i_e(3)^2*q1*q2*q3^2 - 4*Iy*w_i_e(3)^2*q1*q2*q4^2 + 4*Iy*w_i_e(3)^2*q1^2*q3*q4 - 4*Iy*w_i_e(3)^2*q2^2*q3*q4 - 2*Ix*w_i_e(3)*q1*q4*wx - 2*Ix*w_i_e(3)*q2*q3*wx - 2*Ix*w_i_e(3)*q1*q3*wy + 2*Iy*w_i_e(3)*q1*q4*wx + 2*Iy*w_i_e(3)*q2*q3*wx + 2*Ix*w_i_e(3)*q2*q4*wy + 2*Iy*w_i_e(3)*q1*q3*wy - 2*Iy*w_i_e(3)*q2*q4*wy)/Iz;

%NEED vel_x vel_y vel_z
%NEED new r_sc

%ACCELERATION FROM TRANSLATIONAL EOM
acc_x = Fdx/Msc - (GM/norm(r_sc)^3)*x;
acc_y = Fdy/Msc - (GM/norm(r_sc)^3)*y;
acc_z = Fdz/Msc - (GM/norm(r_sc)^3)*z;

%BUILD THE VECTORS
vel   = [vel_x; vel_y; vel_z];
accel = [acc_x; acc_y; acc_z];
omega_dot = [omega_dot_x; omega_dot_y; omega_dot_z];
dq = [dq1; dq2; dq3; dq4];

%RETURN THE DERIVATIVE OF THE STATE
state_dot = [vel; accel; omega_dot; dq];

end %CHOMPTT_EOM.m -------------------------------------------------------
