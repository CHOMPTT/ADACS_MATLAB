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
wa = state(14);% !!!!
wb = state(15);% !!!!
wc = state(16) ;% !!!!

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

%RW INERTIA CONSTANTS
Mw=1;%Reaction Wheel Mass constant in kg %%%REPLACE WITH ACTUAL VALUE
Rw=1;%Reaction Wheel radii constant (m)  %%%REPLACE WITH ACTUAL VALUE

%SPACECRAFT MOMENTS OF INERTIA
Ix=1/12*Msc*(b^2+c^2) ;            
Iy=1/12*Msc*(a^2+c^2) ;
Iz=1/12*Msc*(a^2+b^2) ;

% REACTION WHEEL MOMENT OF INERTIA
Iw = (1/2)*Mw*(Rw^2); % !!!!

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

%ANGULAR ACCELERATION OF SPACECRAFT FROM ROTATIONAL EOM (from SC_RW_EOM2.m)
omega_dot_x = -(Tcx - Tdx - Iw*wb*wz + Iw*wc*wy - Iy*wy*wz + Iz*wy*wz + Iw*We*q1^2*wb + Iw*We*q2^2*wb - Iw*We*q3^2*wb - Iw*We*q4^4*wb + Iy*We*q1^2*wy + Iy*We*q2^2*wy - Iy*We*q3^2*wy - Iy*We*q4^4*wy - Iz*We*q1^2*wy - Iz*We*q2^2*wy + Iz*We*q3^2*wy + Iz*We*q4^4*wy - 2*Iy*We^2*q2*q3^3 + 2*Iy*We^2*q1^3*q4 + 2*Iy*We^2*q2^3*q3 - 2*Iy*We^2*q1*q4^5 + 2*Iz*We^2*q2*q3^3 - 2*Iz*We^2*q1^3*q4 - 2*Iz*We^2*q2^3*q3 + 2*Iz*We^2*q1*q4^5 + 2*Ix*We*dq1*q3 + 2*Ix*We*dq3*q1 - 2*Ix*We*dq2*q4 - 2*Ix*We*dq4*q2 + 2*Iy*We^2*q1^2*q2*q3 + 2*Iy*We^2*q1*q2^2*q4 - 2*Iy*We^2*q1*q3^2*q4 - 2*Iy*We^2*q2*q3*q4^4 - 2*Iz*We^2*q1^2*q2*q3 - 2*Iz*We^2*q1*q2^2*q4 + 2*Iz*We^2*q1*q3^2*q4 + 2*Iz*We^2*q2*q3*q4^4 + 2*Iw*We*q1*q4*wc + 2*Iw*We*q2*q3*wc - 2*Iy*We*q1*q4*wz - 2*Iy*We*q2*q3*wz + 2*Iz*We*q1*q4*wz + 2*Iz*We*q2*q3*wz)/Ix;  % !!!!
omega_dot_y = -(Tcy - Tdy + Iw*wa*wz - Iw*wc*wx + Ix*wx*wz - Iz*wx*wz - Iw*We*q1^2*wa - Iw*We*q2^2*wa + Iw*We*q3^2*wa + Iw*We*q4^4*wa - Ix*We*q1^2*wx - Ix*We*q2^2*wx + Ix*We*q3^2*wx + Ix*We*q4^4*wx + Iz*We*q1^2*wx + Iz*We*q2^2*wx - Iz*We*q3^2*wx - Iz*We*q4^4*wx + 2*Ix*We^2*q1*q3^3 - 2*Ix*We^2*q1^3*q3 + 2*Ix*We^2*q2^3*q4 - 2*Ix*We^2*q2*q4^5 - 2*Iz*We^2*q1*q3^3 + 2*Iz*We^2*q1^3*q3 - 2*Iz*We^2*q2^3*q4 + 2*Iz*We^2*q2*q4^5 + 2*Iy*We*dq1*q4 + 2*Iy*We*dq2*q3 + 2*Iy*We*dq3*q2 + 2*Iy*We*dq4*q1 - 2*Ix*We^2*q1*q2^2*q3 + 2*Ix*We^2*q1^2*q2*q4 - 2*Ix*We^2*q2*q3^2*q4 + 2*Ix*We^2*q1*q3*q4^4 + 2*Iz*We^2*q1*q2^2*q3 - 2*Iz*We^2*q1^2*q2*q4 + 2*Iz*We^2*q2*q3^2*q4 - 2*Iz*We^2*q1*q3*q4^4 - 2*Iw*We*q1*q3*wc + 2*Iw*We*q2*q4*wc + 2*Ix*We*q1*q3*wz - 2*Ix*We*q2*q4*wz - 2*Iz*We*q1*q3*wz + 2*Iz*We*q2*q4*wz)/Iy;  % !!!!
omega_dot_z = (Tdz - Tcz + Iw*wa*wy - Iw*wb*wx + Ix*wx*wy - Iy*wx*wy - 4*Iz*We*dq4*q4^3 + 2*Iz*We*dq1*q1 + 2*Iz*We*dq2*q2 - 2*Iz*We*dq3*q3 + 4*Ix*We^2*q1*q2*q3^2 - 4*Ix*We^2*q1*q2*q4^2 + 4*Ix*We^2*q1^2*q3*q4 - 4*Ix*We^2*q2^2*q3*q4 - 4*Iy*We^2*q1*q2*q3^2 + 4*Iy*We^2*q1*q2*q4^2 - 4*Iy*We^2*q1^2*q3*q4 + 4*Iy*We^2*q2^2*q3*q4 + 2*Iw*We*q1*q4*wa + 2*Iw*We*q2*q3*wa - 2*Iw*We*q1*q3*wb + 2*Iw*We*q2*q4*wb + 2*Ix*We*q1*q4*wx + 2*Ix*We*q2*q3*wx + 2*Ix*We*q1*q3*wy - 2*Iy*We*q1*q4*wx - 2*Iy*We*q2*q3*wx - 2*Ix*We*q2*q4*wy - 2*Iy*We*q1*q3*wy + 2*Iy*We*q2*q4*wy)/Iz;  % !!!!

%ANGULAR ACCELERATION OF REACTION WHEELS FROM ROTATIONAL EOM (from SC_RW_EOM.m)
omega_dot_a = (Tcx - Iw*dwx - 2*Iw*We*dq1*q3 - 2*Iw*We*dq3*q1 + 2*Iw*We*dq2*q4 + 2*Iw*We*dq4*q2)/Iw;  % !!!!
omega_dot_b = -(Iw*dwy - Tcy + 2*Iw*We*dq1*q4 + 2*Iw*We*dq2*q3 + 2*Iw*We*dq3*q2 + 2*Iw*We*dq4*q1)/Iw;  % !!!!
omega_dot_c = (- 4*Iw*We*dq4*q4^3 + Tcz - Iw*dwz + 2*Iw*We*dq1*q1 + 2*Iw*We*dq2*q2 - 2*Iw*We*dq3*q3)/Iw;  % !!!!

%ACCELERATION OF SPACECRAFT FROM TRANSLATIONAL EOM
acc_x = (Fdx/Msc - (GM_E*x)/(x^2 + y^2 + z^2)^(3/2));
acc_y = (Fdy/Msc - (GM_E*y)/(x^2 + y^2 + z^2)^(3/2));
acc_z = (Fdz/Msc - (GM_E*z)/(x^2 + y^2 + z^2)^(3/2));

%BUILD THE VECTORS
vel   = [vel_x; vel_y; vel_z];
accel = [acc_x; acc_y; acc_z];
omega_dot = [omega_dot_x; omega_dot_y; omega_dot_z];
rw_omega_dot = [omega_dot_a; omega_dot_b; omega_dot_c];
dq = [dq1; dq2; dq3; dq4];

%RETURN THE DERIVATIVE OF THE STATE
state_dot = [vel; accel; omega_dot; dq; rw_omega_dot];