% REACTION WHEEL ROTATIONAL DYNAMICS
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
W_i_w = W_i_s + W_s_w;    %omega wheels to inertial

% zero equation for the wheels
% Tc = S_d/dt(H_wheels), where H_wheels = I_Wheels * I_omega_W
H_wheel = Jw*W_i_w;
H_satellite = J*W_i_s;
z1 = -Tc + diff(H_wheel,t);

% convert time to symbolic so that we can solve for angular acceleration of
% the weels, dwa, dwb, dwc explicitly
clear wx wy wz wa wb wc q1 q2 q3 q4 % clear the time-variables
syms wx wy wz wa wb wc q1 q2 q3 q4  % declare the variables as symbols
syms dwx dwy dwz dwa dwb dwc dq1 dq2 dq3 dq4
z1_sym(1,1) = eval(RW_convTimeToSym(z1(1)));
z1_sym(2,1) = eval(RW_convTimeToSym(z1(2)));
z1_sym(3,1) = eval(RW_convTimeToSym(z1(3)));

% solve z1 for dwa dwb dwc
dW_wheel = solve(z1_sym(1)==0,z1_sym(2)==0,z1_sym(3)==0,dwa,dwb,dwc )
% !! your solution will be in
% !! dW_wheel.dwa, dW_wheel.dwb, dW_wheel.dwc
% !! copy and paste the solutions into SC_RW_EOM2.m for the spacecraft dynamics
% !! and CHOMPTT_EOM.m
dW_wheel.dwa
dW_wheel.dwb
dW_wheel.dwc
