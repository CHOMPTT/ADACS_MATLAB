% Plot_State.m -----------------------------------------------------------
% 
% This file is used to plot the state vector
% INPUT:
%   state is the state vector (13 x n)
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
%   t = 
%
% OUTPUT:
%   None
%
%%%Plotting state vectors versus time
%%%Each state variable has its own window  
function Plot_State(state_vec, t);
   for k=1:13
	figure(k)
	plot(t(1:numel(t)-1),state_vec(k,:));
	xlabel('Time (s)')
	if k==1||k==2||k==3
    	ylabel('Position(m)')
    	if k == 1
        	title('X Position')
    	elseif k==2
        	title('Y Position')
    	elseif k==3
        	title ('Z Position')
    	end
	elseif k==4|| k==5||k==6
    	ylabel('Velocity(m/s)')
    	if k==4
        	title('X Velocity')
    	elseif k==5
        	title('Y Velocity')
    	elseif k==6
        	title('Z Velocity')
    	end
	elseif k==7||k==8||k==9
    	ylabel('Angular Velocity(rad/s)')
    	if k==7
        	title('X Angular Velocity')
    	elseif k==8
        	title('Y Angular Velocity')
    	elseif k==9
        	title('Z Angular Velocity')
    	end
	elseif k==10||k==11||k==12||k==13
    	ylabel('Quaternion')
    	if k==10
        	title('q1')
    	elseif k==11
        	title('q2')
    	elseif k==12
        	title('q3')
    	elseif k==13
        	title('q4')
    	end
	end
   end
end