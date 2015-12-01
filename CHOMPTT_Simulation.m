%main.m -------------------------------------------------------------------
clear all;

%VARIABLE(S) INITIALIZATION
%Prev_B_Vect = [0, 0, 0];
wE = 0.7292115 * 10^-4;  %Earth's angular velocity in (rad/s) ECI/ECEF  
GM_E = 3.986004415e+14;   % (m3/s2)

w_e_s = [0 ; 0; 0];                 %S/C angular velocity w.r.t ECEF in S/C frame assumed given by sensors
q = [0; 0; 0; 1];                   %Quaternion matrix, assumed given by sensors

R= 6378137;                 %radius of the earth, m
height = 400*1000;          % altitude of satellite, m
r_initial = R + height;     % Initial radius of the satellite in orbit (m)
r_sc = [0; 0; r_initial];   % position of satellite in ECI  
vCirc = sqrt(GM_E/r_initial); % Get initial velocity of a circular orbit
vel_sc = [vCirc; 0; 0];      % tangential velocity of satellite in ECI
T = 2*pi*sqrt(((r_initial)^3)/GM_E); % sec << Calculate orbit period
    
%TIME VECTOR (s)
t_start = 1;
t_end = 3*T;    % 3 orbits
t_step = 1;
t = t_start:t_step:t_end;%time vector in seconds

Fd_vec = importdata('ElaNaXIX_DisturbanceForces.mat');% Disturbance force vector in ECI ref. frame, in Newtons

%Sat inertia constants
Msc=3;      %Spacecraft Mass := 3Kg
a=0.1;      %Spacecraft length along x-axis  (m)
b=0.1;      %Spacecraft length along y-axis  (m)
c=0.34;     %Spacecraft length along z-axis  (m)
% center of mass 
Rx = a/2 ;
Ry = b/2 ;
Rz = c/2 ;
Rcom_i = [Rx ; Ry; Rz];% distance from c.o.m to side of s/c in direction of unit vector 
        
%State Vector
state = [r_sc ; vel_sc ; w_e_s ; q];

%state_vec = zeros(1:length(t));
%state(1)=r_sc(1)=x position of satellite in ECI
%state(2)=r_sc(2)=y position of satellite in ECI
%state(3)=r_sc(3)=z position of satellite in ECI
%state(4)=vel_sc(1)= x (along-track) velocity of satellite in ECI
%state(5)=vel_sc(2)= y (cross-track) velocity of satellite in ECI
%state(6)=vel_sc(3)= z (radial) velocity of satellite in ECI
%state(7)=w_e_s(1)= x S/C angular velocity w.r.t ECEF in S/C frame
%state(8)=w_e_s(2)= x S/C angular velocity w.r.t ECEF in S/C frame
%state(9)=w_e_s(3)= x S/C angular velocity w.r.t ECEF in S/C frame
%state(10)=q(1)
%state(11)=q(2)
%state(12)=q(3)
%state(13)=q(4)


%BEGIN SIMULATION
for i=1:length(t) 

        %ROTATION MATRICES 
        %ECI -> ECEF : Inertial to Earth-Fixed Ref. Frame Rotation
        theta_earth = wE*t(i);    % angle of rotation (rad)
        R_e_i = [ cos(theta_earth), sin(theta_earth), 0; 
                 -sin(theta_earth), cos(theta_earth), 0; 
                  0,                0,                1];
        %ECEF -> S/C : Earth Fixed to Spacecraft Fixed Ref. Frame
        q = [state(10); state(11); state(12); state(13)];   % Set the current quaternion
        R_s_e = [q(1)^2-q(2)^2-q(3)^2+q(4)^2, 2*(q(1)*q(2)+q(3)*q(4)), 2*(q(1)*q(3)-q(2)*q(4)); 
                 2*(q(1)*q(2)-q(3)*q(4)), -q(1)^2+q(2)^2-q(3)^2+q(4)^2, 2*(q(2)*q(3)+q(1)*q(4)); 
                 2*(q(1)*q(3)+q(2)*q(4)), 2*(q(2)*q(3)+q(1)*q(4)), -q(1)^2-q(2)^2+q(3)^2+q(4)^4] ;
        %ECI -> S/C : Inertial to Spacecraft (Body) Fixed Ref. Frame
        R_s_i = R_s_e*R_e_i;

        %DISTURBANCE FORCE & TORQUE BLOCK 
%         Fd_s = R_s_i * Fd_vec(:,i);     % Disturbance force expressed in S/C ref. frame
%         Td = cross(Rcom_i, Fd_s);   % Compute the disturbance Torque from Fd expressed in S/C Ref. Frame
        
        %EPS/MODE BLOCK 
        %Select a mode for the system you want to run (RW or MT) by uncommenting
          %EPS_MODE = 'POINTING';
          %EPS_MODE = 'DETUMBLE';
          EPS_MODE = 'DRIFTING';

        %CONTROLS BLOCK 
         if EPS_MODE == 'DETUMBLE';
             k = 4*10^4; %Proportional Gain
             MAX_DIPOLE = 0.03774; %Am^2 %Hardware Limitation
             [Tc , Prev_B_Vect] = B_DOT_CONTROL(k, state, MAX_DIPOLE, i, t_step, Prev_B_Vect);
         end

         if EPS_MODE == 'POINTING';
            k = [0.2;0.2 ;0.2];%Proportional Gain
            c = [0.5;0.5; 0.5];%Derivative Gain
            MAX_TORQUE = 3.75e-3; %N*m %Hardware Limitation
            Tc = RW_PD_CONTROL(k, c, state, MAX_TORQUE);
         end

         if EPS_MODE == 'DRIFTING';
            Tc = [0;0;0];
         end

        %ACTUATION NOISE BLOCK 
        %Add actuation noise to the commanded torque as an 'actual torque'
         if EPS_MODE == 'DETUMBLE';
            sigma = 0;
         end

         if EPS_MODE == 'POINTING';
            sigma = 0;
         end

         if EPS_MODE == 'DRIFTING';
            sigma = 0;
         end
         Tc = normrnd(Tc,sigma,3,1);

        %MEASUREMENT NOISE BLOCK 
        %Add sensor noise to the 'true' state from the dynamics block 

        %DYNAMICS BLOCK
        %RK4 Integrator
        Td = zeros(3,1);
        Tc = zeros(3,1);
        Fd = zeros(3,1);
            k1 = CHOMPTT_EOM(t(i)           , state              , Tc, Fd, Td);     
            k2 = CHOMPTT_EOM(t(i) + t_step/2, state + t_step*k1/2, Tc, Fd, Td); 
            k3 = CHOMPTT_EOM(t(i) + t_step/2, state + t_step*k2/2, Tc, Fd, Td); 
            k4 = CHOMPTT_EOM(t(i)           , state + t_step*k3  , Tc, Fd, Td); 
            state = state + (1/6)*t_step*(k1+2*k2+2*k3+k4); % Generate the state at the current time-step
            state_vec(:,i) = state; % Save the state vector for later
end % end of simulation

figure, plot(t./3600, state_vec(1:3,:)./1000), grid on,...
        title('SC position'),...
        xlabel('time (hr)'),...
        ylabel('position (km)'),...
        legend('x','y','z');

%end of main.m ------------------------------------------------------------------
