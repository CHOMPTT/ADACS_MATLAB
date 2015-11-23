%main.m -------------------------------------------------------------------
clear all;

%TIME VECTOR (s)
t_start = 1;
t_end = 10001;
t_end = 5400*5;
t_step = 10;
t = t_start:t_step:t_end;%time vector in seconds

%VARIABLE(S) INITIALIZATION
%Prev_B_Vect = [0, 0, 0];
w_i_e = [0; 0; 0.7292115 * 10^-4];%Earth's angular velocity in ECI/ECEF  
w_e_s = [0 ; 0; 0];%S/C angular velocity w.r.t ECEF in S/C frame assumed given by sensors
q = [0; 0; 0; 1];%Quaternion matrix, assumed given by sensors

R= 6378137;%radius of the earth, m
height = 400*1000;% altitude of satellite, m
r_initial = R + height;% Initial radius of the satellite in orbit (m)
r_sc = [0; 0; r_initial];% position of satellite in ECI  
vel_sc = [0; 0 ; 0];% tangential velocity of satellite in ECI

%Sat inertia constants
Msc=3;%Spacecraft Mass := 3Kg
a=0.1;%Spacecraft length along x-axis  (m)
b=0.1;%Spacecraft length along y-axis  (m)
c=0.34;%Spacecraft length along z-axis  (m)

%State Vector
state = [r_sc ; vel_sc ; w_e_s; q];
state_vec = zeros(1:length(t));
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


%BEGIN SIMULATION
for i=1:length(t)-1  
tsim = t(i);
    for j = 1:4
        if j == 1
            tstate = tsim;
        elseif j == 2
            tstate = tsim+t_step/2;
        elseif j == 3
            tstate = tsim+t_step/2;
        elseif j == 4
            tstate = tsim+t_step;
        end

        %ROTATION MATRICES 
        %ECI -> ECEF : Inertial to Earth-Fixed Ref. Frame Rotation
        theta_earth = w_i_e(3)*tstate;% angle of rotation (rad)
        R_e_i = [cos(theta_earth), sin(theta_earth), 0; -sin(theta_earth), cos(theta_earth), 0; 0, 0, 1];
        %ECEF -> S/C : Earth Fixed to Spacecraft Fixed Ref. Frame 
        R_s_e = [q(1)^2-q(2)^2-q(3)^2+q(4)^2, 2*(q(1)*q(2)+q(3)*q(4)), 2*(q(1)*q(3)-q(2)*q(4)); 
                 2*(q(1)*q(2)-q(3)*q(4)), -q(1)^2+q(2)^2-q(3)^2+q(4)^2, 2*(q(2)*q(3)+q(1)*q(4)); 
                 2*(q(1)*q(3)+q(2)*q(4)), 2*(q(2)*q(3)+q(1)*q(4)), -q(1)^2-q(2)^2+q(3)^2+q(4)^4] ;
        %ECI -> S/C : Inertial to Spacecraft (Body) Fixed Ref. Frame
        R_s_i = R_s_e*R_e_i;

        %DISTURBANCE FORCE & TORQUE BLOCK 
        %Get the disturbance at the current time step given the disturbance
        %force vector Fd (3x1)                                
        %c.o.m = center of mass
        c = 0.034;
        Rx = a/2 ;
        Ry = b/2 ;
        Rz = c/2 ;
        Rcom_i = [Rx ; Ry; Rz];% distance from c.o.m to side of s/c in direction of unit vector 
        Fd = importdata('ElaNaXIX_DisturbanceForces.mat');% Disturbance force vector in ECI ref. frame, in Newtons
        Fd_s = R_s_i* Fd(:,i);%Disturbance force expressed in S/C ref. frame
        Td = cross(Rcom_i, Fd_s);%Compute the disturbance Torque from Fd expressed in S/C Ref. Frame

        %EPS/MODE BLOCK 
        %Select a mode for the system you want to run (RW or MT) by uncommenting
          %EPS_MODE = 'POINTING';
          %EPS_MODE = 'DETUMBLE';
          EPS_MODE = 'DRIFTING';

        %CONTROLS BLOCK 
         if EPS_MODE == 'DETUMBLE';
             k = 4*10^4;%Proportional Gain
             MAX_DIPOLE = 0.03774; %Am^2 %Hardware Limitation
             [Tc , Prev_B_Vect] = B_DOT_CONTROL(k, state, MAX_DIPOLE, i, t_step, Prev_B_Vect);
         end

         if EPS_MODE == 'POINTING';
            k= [0.2;0.2 ;0.2];%Proportional Gain
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
        if j == 1
            k1 = CHOMPTT_EOM(state, Tc, Fd, Td, w_i_e);     
        elseif j == 2
            k2 = CHOMPTT_EOM(state + t_step*k1/2, Tc, Fd, Td, w_i_e); 
        elseif j == 3
            k3 = CHOMPTT_EOM(state + t_step*k2/2, Tc, Fd, Td, w_i_e); 
        elseif j == 4
            k4 = CHOMPTT_EOM(state + t_step*k3  , Tc, Fd, Td, w_i_e);  
        end 
    end 

    state = state + (1/6)*t_step*(k1+2*k2+2*k3+k4);%Generate the state at the current time-step
    state_vec(:,i) = state;%Save the state vector for later
end % end of simulation

%Plot your results

%end of main.m ---------------------------------------------------
