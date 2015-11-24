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

function [state_dot] = CHOMPTT_EOM(t, state, Tc, Fd, Td)     

%GET STATE VARIABLES
x  =  state(1);        
y  =  state(2);        
z  =  state(3);         
vel_x =  state(4);         
vel_y =  state(5);         
vel_z =  state(6);        
wx =  state(7);        
wy =  state(8);         
wz =  state(9);         
q1 = state(10);        
q2 = state(11);        
q3 = state(12);         
q4 = state(13);   

% Postion Vector in the inertial reference frame
r_sc = state(1:3);
% Quaternion vector ijk
q = [q1;q2;q3];
% Angular velocity vector
omega = [wx; wy; wz];
% Angular velocity of the Earth
wE = 0.7292115e-4; % rad/sec
%CONSTANT(s)
GM_E = 3.986004415e+14;   % (m3/s2)

%SAT INERTIA CONSTANTS
Msc=3;%Spacecraft Mass := 3Kg
a=0.1;%Spacecraft length along x-axis  (m)
b=0.1;%Spacecraft length along y-axis  (m)
c=0.34;%Spacecraft length along z-axis  (m)

%SPACECRAFT MOMENTS OF INERTIA
Ix=1/12*Msc*(b^2+c^2) ;            
Iy=1/12*Msc*(a^2+c^2) ;
Iz=1/12*Msc*(a^2+b^2) ;

%CONTROL FORCES
Tcx = Tc(1);
Tcy = Tc(2);
Tcz = Tc(3);

%DISTURBANCE FORCES AND TORQUES
Fdx = Fd(1); % Disturbance forces in the Inertial reference frame
Fdy = Fd(2);    
Fdz = Fd(3); 
Tdx = Td(1); % Disturbance torques in the spacecraft reference frame
Tdy = Td(2);
Tdz = Td(3);

%QUATERNION DERIVATIVES (EOM GIVEN)
dq = 0.5*cross(omega,q) + 0.5*q4*omega;
dq1 = dq(1) ;
dq2 = dq(2);
dq3 = dq(3);
dq4 = -0.5*dot(omega,q);

%ANGULAR ACCELERATION FROM ROTATIONAL EOM
omega_dot_x = (Tcx + Tdx + Iy*wy*wz - Iz*wy*wz - Iy*q1^2*wE*wy - Iy*q2^2*wE*wy + Iy*q3^2*wE*wy + Iy*q4^2*wE*wy + Iz*q1^2*wE*wy + Iz*q2^2*wE*wy - Iz*q3^2*wE*wy - Iz*q4^2*wE*wy + 2*Iy*q1*q4^3*wE^2 + 2*Iy*q2*q3^3*wE^2 - 2*Iy*q1^3*q4*wE^2 - 2*Iy*q2^3*q3*wE^2 - 2*Iz*q1*q4^3*wE^2 - 2*Iz*q2*q3^3*wE^2 + 2*Iz*q1^3*q4*wE^2 + 2*Iz*q2^3*q3*wE^2 - 2*Ix*dq1*q3*wE - 2*Ix*dq3*q1*wE + 2*Ix*dq2*q4*wE + 2*Ix*dq4*q2*wE - 2*Iy*q1^2*q2*q3*wE^2 - 2*Iy*q1*q2^2*q4*wE^2 + 2*Iy*q1*q3^2*q4*wE^2 + 2*Iy*q2*q3*q4^2*wE^2 + 2*Iz*q1^2*q2*q3*wE^2 + 2*Iz*q1*q2^2*q4*wE^2 - 2*Iz*q1*q3^2*q4*wE^2 - 2*Iz*q2*q3*q4^2*wE^2 + 2*Iy*q1*q4*wE*wz + 2*Iy*q2*q3*wE*wz - 2*Iz*q1*q4*wE*wz - 2*Iz*q2*q3*wE*wz)/Ix;
omega_dot_y = -(Ix*wx*wz - Tdy - Tcy - Iz*wx*wz - Ix*q1^2*wE*wx - Ix*q2^2*wE*wx + Ix*q3^2*wE*wx + Ix*q4^2*wE*wx + Iz*q1^2*wE*wx + Iz*q2^2*wE*wx - Iz*q3^2*wE*wx - Iz*q4^2*wE*wx + 2*Ix*q1*q3^3*wE^2 - 2*Ix*q1^3*q3*wE^2 - 2*Ix*q2*q4^3*wE^2 + 2*Ix*q2^3*q4*wE^2 - 2*Iz*q1*q3^3*wE^2 + 2*Iz*q1^3*q3*wE^2 + 2*Iz*q2*q4^3*wE^2 - 2*Iz*q2^3*q4*wE^2 + 2*Iy*dq1*q4*wE + 2*Iy*dq2*q3*wE + 2*Iy*dq3*q2*wE + 2*Iy*dq4*q1*wE - 2*Ix*q1*q2^2*q3*wE^2 + 2*Ix*q1^2*q2*q4*wE^2 + 2*Ix*q1*q3*q4^2*wE^2 - 2*Ix*q2*q3^2*q4*wE^2 + 2*Iz*q1*q2^2*q3*wE^2 - 2*Iz*q1^2*q2*q4*wE^2 - 2*Iz*q1*q3*q4^2*wE^2 + 2*Iz*q2*q3^2*q4*wE^2 + 2*Ix*q1*q3*wE*wz - 2*Ix*q2*q4*wE*wz - 2*Iz*q1*q3*wE*wz + 2*Iz*q2*q4*wE*wz)/Iy;
omega_dot_z = (Tcz + Tdz + Ix*wx*wy - Iy*wx*wy + 2*Iz*dq1*q1*wE + 2*Iz*dq2*q2*wE - 2*Iz*dq3*q3*wE - 2*Iz*dq4*q4*wE + 4*Ix*q1*q2*q3^2*wE^2 - 4*Ix*q1*q2*q4^2*wE^2 + 4*Ix*q1^2*q3*q4*wE^2 - 4*Ix*q2^2*q3*q4*wE^2 - 4*Iy*q1*q2*q3^2*wE^2 + 4*Iy*q1*q2*q4^2*wE^2 - 4*Iy*q1^2*q3*q4*wE^2 + 4*Iy*q2^2*q3*q4*wE^2 + 2*Ix*q1*q4*wE*wx + 2*Ix*q2*q3*wE*wx + 2*Ix*q1*q3*wE*wy - 2*Iy*q1*q4*wE*wx - 2*Iy*q2*q3*wE*wx - 2*Ix*q2*q4*wE*wy - 2*Iy*q1*q3*wE*wy + 2*Iy*q2*q4*wE*wy)/Iz;

%ACCELERATION FROM TRANSLATIONAL EOM
acc_x = (Fdx - (GM_E/((norm(r_sc))^3))*x)/Msc;
acc_y = (Fdy - (GM_E/((norm(r_sc))^3))*y)/Msc;
acc_z = (Fdz - (GM_E/((norm(r_sc))^3))*z)/Msc;


%BUILD THE VECTORS
vel   = [vel_x; vel_y; vel_z];
accel = [acc_x; acc_y; acc_z];
omega_dot = [omega_dot_x; omega_dot_y; omega_dot_z];
dq = [dq1; dq2; dq3; dq4];

%RETURN THE DERIVATIVE OF THE STATE
state_dot = [vel; accel; omega_dot; dq];


