 % main.m ------------------------------------------------------------------
clear all;
% state your time vector  (s)

t_start = 1;
t_end =10001;
t_end =5400*5
t_step = 10;

t=t_start:t_step:t_end ;    	% time vector in seconds

% initialize your variables
B_Vect = 1.0e+04*[0, 0, 0];
B_Vect = transpose(B_Vect);

w_i_e = [0; 0; 0.7292115 * 10^-4] ;   % Earth’s angular velocity in ECI/ECEF  
w_e_s = [0 ; 0; 0] ;		           % S/C angular velocity w.r.t ECEF in S/C frame
				           % assumed given by sensors
q = [0; 0; 0; 1] ;                                  % Quaternion matrix, assumed given by sensors


R= 6378137;                            	 % radius of the earth, m
height = 400*1000;                    	 % altitude of satellite, m
r_initial = R + height;                           % Initial radius of the satellite in orbit (m)
r_sc = [0; 0; r_initial] ;       	             % position of satellite in ECI  
vel_sc = [0; 0 ; 0] ;			 % tangential velocity of satellite in ECI

% Sat inertia constants
Msc=3;     %Spacecraft Mass := 3Kg
a=0.1;      %Spacecraft length along x-axis  (m)
b=0.1;      %Spacecraft length along y-axis  (m)
c=0.34;    %Spacecraft length along z-axis  (m)

 state = [r_sc ; vel_sc ; w_e_s; q] ;
%state(1)=r_sc(1)=x position of satellite in ECI
%state(2)=r_sc(2)=y position of satellite in ECI
%state(3)=r_sc(3)=z position of satellite in ECI
%state(4)=vel_sc(1)= x tangential velocity of satellite in ECI
%state(5)=vel_sc(2)= y tangential velocity of satellite in ECI
%state(6)=vel_sc(3)= z tangential velocity of satellite in ECI
%state(7)=w_e_s(1)= x S/C angular velocity w.r.t ECEF in S/C frame
%state(8)=w_e_s(2)= x S/C angular velocity w.r.t ECEF in S/C frame
%state(9)=w_e_s(3)= x S/C angular velocity w.r.t ECEF in S/C frame
%state(10)=q(1)
%state(11)=q(2)
%state(12)=q(3)
%state(13)=q(4)
% begin simulation
    for i=1:length(t)-1  
       
% ROTATION MATRICES 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% ECI -> ECEF : Inertial to Earth-Fixed Ref. Frame Rotation %%%%%%

theta_earth = w_i_e(3)*t(i);                               % angle of rotation (rad)
R_e_i = [cos(theta_earth), sin(theta_earth), 0; -sin(theta_earth), cos(theta_earth), 0; 0, 0, 1];


%%%%%%% ECEF -> S/C : Earth Fixed to Spacecraft Fixed Ref. Frame %%%%%

R_s_e = [q(1)^2-q(2)^2-q(3)^2+q(4)^2, 2*(q(1)*q(2)+q(3)*q(4)), 2*(q(1)*q(3)-q(2)*q(4)); 2*(q(1)*q(2)-q(3)*q(4)), -q(1)^2+q(2)^2-q(3)^2+q(4)^2, 2*(q(2)*q(3)+q(1)*q(4)); 2*(q(1)*q(3)+q(2)*q(4)), 2*(q(2)*q(3)+q(1)*q(4)), -q(1)^2-q(2)^2+q(3)^2+q(4)^4] ;

%%%%% ECI -> S/C : Inertial to Spacecraft (Body) Fixed Ref. Frame %%%%
R_s_i = R_s_e*R_e_i   ;


% DISTURBANCE FORCE & TORQUE BLOCK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Get the disturbance at the current time step given the disturbance
    %   force vector Fd (3x1)                                 %c.o.m = center of max
    %
     c=0.034;
     Rx = a/2 ;
     Ry = b/2 ;
     Rz = c/2 ;
     
     Rcom_i = [Rx ; Ry; Rz];             % distance from c.o.m to side of s/c in direction of unit vector 
     
     Fd = importdata('ElaNaXIX_DisturbanceForces.mat');  % Disturbance force vector in ECI ref. frame, in Newtons
                
   % Disturbance force expressed in S/C ref. frame
     Fd_s = R_s_e *R_e_i* Fd(:,i)  ;

   %  Compute the disturbance Torque from Fd
    Td = cross(Rcom_i, Fd_s) ;                  % expressed in S/C Ref. Frame
    %Td = [0;0;0];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % EPS/MODE BLOCK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %   Select a mode for the whichever system you want to run (RW or MT)
EPS_MODE='POINTING';%DETUMBLE or POINTING
%EPS_MODE='DETUMBLE';%DETUMBLE or POINTING

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
 % CONTROLS BLOCK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Get the current Euler Angles from the quaternion
    qt = transpose(state(10:13));
    zyx = Quat2Eul321(qt);
    
    %   Apply your control law (PD and B-dot) here
    %   Compute commanded torque 
    %
if EPS_MODE == 'DETUMBLE';
% %BDOT
% %% convert the current R-vector from S/C into ECEF from there use ECEF Position to LLA 
r_sct = transpose(state(1:3));
 lla = ecef2lla(r_sct) ; 
 
 [current_B_Vect, H, DEC, DIP, F] = wrldmagm(lla(3), lla(1), lla(2), decyear(2015,10,25),'2015');
 %disp(current_B_Vect)
 B_Vect(:,i) = [current_B_Vect*10^-09];
 %B_Vect = transpose(B_Vect)
     K = 4*10^4; %Proportional Gain
     %K = 1;
     if (i<2);
     B_DOT(:,1)=[0,0,0];
     else
         B_DOT(:,i)= (B_Vect(:,i)-B_Vect(:,i-1))*(1/t_step);
     end
     M = -K*B_DOT(:,i);
     if M == 0
     M = [0,0,0];
     end
     %limit for max dipole moment
     if size(M)==[1,3]
     M=transpose(M);
     end
     max_dipole_moment = 0.03774; %Am^2
     for j = 1:length(M);
         if (M(j,1)) > max_dipole_moment
 	       M(j,1) = max_dipole_moment;
         end
        if (M(j,1)) < -max_dipole_moment
 	       M(j,1) = -max_dipole_moment;
        end
     end
 disp(M)
 B_Vect(:,i)= transpose(B_Vect(:,i));
 Tc=  cross(M,B_Vect(:,i));
end





if EPS_MODE == 'POINTING'

%PD Controller for Reaction Wheels
k = [0.2;0.2 ;0.2 ]; 
c = [0.5;0.5; 0.5 ];
max_torque = 3.75e-3;

T_x = -k(1,1)*(zyx(3)) - c(1,1)*(state(7));
T_y = -k(2,1)*(zyx(2)) - c(2,1)*(state(8));
T_z = -k(3,1)*(3.1415-abs(zyx(1))) - c(3,1)*(state(8));
T = [T_x; T_y; T_z];

% Checking for Torque Max-Min Requirements
count = length(T);
for j = 1:count
   if abs(T(j,1)) >= max_torque
        if sign(T(j,1)) == 1
            T(j,1) = max_torque;
        else
            T(j,1) = -max_torque;
        end
    end
end
Tc = T;
disp(Tc-Td)
end

%Tc = [0;0;0];
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % ACTUATION NOISE BLOCK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %   Add actuation noise to the commanded torque as an 'actual torque'
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % DYNAMICS BLOCK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    % RK4 Integrator
    k1 = CHOMPTT_EOM(         t, state,               Tc, Fd, Td, w_i_e);         
    k2 = CHOMPTT_EOM(t+t_step/2, state + t_step*k1/2, Tc, Fd, Td, w_i_e); 
    k3 = CHOMPTT_EOM(t+t_step/2, state + t_step*k2/2, Tc, Fd, Td, w_i_e); 
    k4 = CHOMPTT_EOM(  t+t_step, state + t_step*k3,   Tc, Fd, Td, w_i_e); 

    % Generate the state at the current time-step
    state = state + (1/6)*t_step*(k1+2*k2+2*k3+k4);
    %disp(state)
    % Save the state vector for later
    state_vec(:,i) = state; 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % MEASUREMENT NOISE BLOCK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %  Add sensor noise to the 'true' state from the dynamics block 
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end % end of simulation

% Plot your results

% end of main.m ------------------------------------------------------------------