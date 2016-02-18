% SPACECRAFT ROTATIONAL DYNAMICS w/ RW
% A.Nguyen 1/25/2016

clear all
syms t

% Angular velocity of the Spacecraft as a function of time
wx = sym('wx(t)')
wy = sym('wy(t)')
wz = sym('wz(t)')
W_e_s = [wx; wy; wz]

% Angular velocity of the individual wheels as a function of time
wa = sym('wa(t)')
wb = sym('wb(t)')
wc = sym('wc(t)')
W_s_w = [wa; wb; wc]

% Disturbance Torques = r_com x Fd
syms Tdx Tdy Tdz
Td = [Tdx; Tdy; Tdz]
% Commanded Control Torques from PID controller
syms Tcx Tcy Tcz
Tc = [Tcx; Tcy; Tcz]
% Inertia matrix of the spacecraft
syms Ix Iy Iz 
J = diag([Ix, Iy, Iz])
% Inertia matrix of the wheels
syms Iw
Jw = diag([Iw, Iw, Iw])   %Moment of Inertia of wheels

% Quaternion of the spacecraft that describes the orientation of the
% spacecraft with respect to the Earth-fixed axis
q1 = sym('q1(t)')
q2 = sym('q2(t)')
q3 = sym('q3(t)')
q4 = sym('q4(t)')
q = [q1, q2, q3, q4]

% Rotation matrix from the Earth fixed to Spacecraft ref. frame
R_s_e = [q(1)^2-q(2)^2-q(3)^2+q(4)^2, 2*(q(1)*q(2)+q(3)*q(4)), 2*(q(1)*q(3)-q(2)*q(4)); 
         2*(q(1)*q(2)-q(3)*q(4)), -q(1)^2+q(2)^2-q(3)^2+q(4)^2, 2*(q(2)*q(3)+q(1)*q(4)); 
         2*(q(1)*q(3)+q(2)*q(4)), 2*(q(2)*q(3)+q(1)*q(4)), -q(1)^2-q(2)^2+q(3)^2+q(4)^4] ;

% Angular velocity of the Earth
syms We
W_i_e = [0;0;We];      
% Angular velocity of the Spacecraft w.r.t. the inertial ref. frame {I}
W_i_s = R_s_e*W_i_e + W_e_s;     
% Angular velocity of the wheels w.r.t. the spacecraft-body ref. frame {S}
W_i_w = W_i_s + W_s_w;    

% zero equation for the wheels
% Tc = S_d/dt(H_wheels), where H_wheels = I_Wheels * I_omega_W
H_wheel = Jw*W_i_w;
H_satellite = J*W_i_s;

% Get the derivatives of the angular momentums by applying transport
% theorem 
dH_wheel = diff(H_wheel,t) + cross(W_i_s, H_wheel)
dH_satellite = diff(H_satellite) + cross(W_i_s, H_satellite)
z2 = -Td + dH_satellite + dH_wheel

% convert time to symbolic so that we can solve for angular acceleration of
% the weels, dwa, dwb, dwc explicitly
clear wx wy wz wa wb wc q1 q2 q3 q4 % clear the time-variables
syms wx wy wz wa wb wc q1 q2 q3 q4  % delcare the variables as symbols
syms dwx dwy dwz dwa dwb dwc dq1 dq2 dq3 dq4
z2_sym(1,1) = eval(RW_convTimeToSym(z2(1)));
z2_sym(2,1) = eval(RW_convTimeToSym(z2(2)));
z2_sym(3,1) = eval(RW_convTimeToSym(z2(3)));

% Substitue the values fo' dwa, dwb, and dwc from SC_RW_EOM.m
clear dwa dwb dwc
% paste the results from SC_RW_EOM.m for dwa, dwb, and dwc here
dwa = (Tcx - Iw*dwx - 2*Iw*We*dq1*q3 - 2*Iw*We*dq3*q1 + 2*Iw*We*dq2*q4 + 2*Iw*We*dq4*q2)/Iw;   % !!!!
dwb = -(Iw*dwy - Tcy + 2*Iw*We*dq1*q4 + 2*Iw*We*dq2*q3 + 2*Iw*We*dq3*q2 + 2*Iw*We*dq4*q1)/Iw;   % !!!!
dwc = (- 4*Iw*We*dq4*q4^3 + Tcz - Iw*dwz + 2*Iw*We*dq1*q1 + 2*Iw*We*dq2*q2 - 2*Iw*We*dq3*q3)/Iw;   % !!!!

z2_sym2(1,1) = eval(z2_sym(1,1));
z2_sym2(2,1) = eval(z2_sym(2,1));
z2_sym2(3,1) = eval(z2_sym(3,1));

% solve z1 for dwa dwb dwc
dW_satellite = solve(z2_sym2(1)==0,z2_sym2(2)==0,z2_sym2(3)==0,dwx,dwy,dwz )
% !! your solution will be in
% !! dW_satellite.dwx, dW_satellite.dwy, dW_satellite.dwz
% !! copy and paste the following ouptut into CHOMPTT_EOM.m
dW_satellite.dwx
dW_satellite.dwy
dW_satellite.dwz



