%%%%%%%%%%%%%%%%% Quaternion to Euler  %%%%%%%%%%%%%%%%%
function zyx = Quat2Eul321(q)
%% Convert quaternions into 321 euler sequence 
% [n,m] = size(q);
% if n > mtest
%     l = n;
% else
%     l = m;
% end
% ROLL = zeros(l,1);
% PITCH = zeros(l,1);
% YAW = zeros(l,1);

% for i = 1:l
    A = 2.*((q(1,3).*q(1,4)) + (q(1,1).*q(1,2)));
    B = 1/2-(q(1,2).^2 + q(1,3).^2);
    ROLL = atan2(A,B);
    
    C = -2.*(-q(1,2).*q(1,4) - q(1,1).*q(1,3));
    PITCH = asin(C);
    
    D = 2*((q(1,2).*q(1,3)) - (q(1,1)*q(1,4)));
    E = 1/2 - (q(1,3).^2 + q(1,4)^2);
    YAW = atan2(D,E);
    
  zyx=[YAW,PITCH,ROLL];
% end

% ROLL = radtodeg(ROLL)
% PITCH = radtodeg(PITCH)
% YAW = radtodeg(YAW)
end
