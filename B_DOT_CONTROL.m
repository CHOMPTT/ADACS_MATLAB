%B_DOT_CONTROL------------------------------------------
% 
% This file is used to simulate the magnatorquer detumble controller (B dot)
%
% INPUT:
%   k is the proportional gain (scalar) 
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
%   MAX_DIPOLE = maximum dipole induced by magantorquer
%
% OUTPUT:
%   Tc is the control torque (3 x 1)
% 

% %% convert the current R-vector from S/C into ECEF from there use ECEF Position to LLA 
function Tc = RW_PD_CONTROL(k, state, MAX_DIPOLE);
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
     %MAX_DIPOLE = 0.03774; %Am^2
     for j = 1:length(M);
         if (M(j,1)) > MAX_DIPOLE
 	       M(j,1) = MAX_DIPOLE;
         end
        if (M(j,1)) < -MAX_DIPOLE
 	       M(j,1) = -MAX_DIPOLE;
        end
     end
 disp(M)
 B_Vect(:,i)= transpose(B_Vect(:,i));
 Tc=  cross(M,B_Vect(:,i));
end





