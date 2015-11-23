function [zyx] = Quat2Eul321(q)
%% Convert quaternions into 321 euler sequence 
%https://en.wikiversity.org/wiki/PlanetPhysics/Direction_Cosine_Matrix_to_E
%uler_321_Angles
%
% %Solving for a R_321
% %Used to solve for the Euler angles in the Quat2Eu321.m
% syms theta phi psi
% %AboutX
% R_1 = [1,0,0;
%        0,cos(psi),sin(psi);
%        0,-sin(psi),cos(psi)];
% %AboutY
% R_2 = [cos(theta),0,-sin(theta);
%       0,1,0;
%       sin(theta),0,cos(theta)];
% %AboutZ
% R_3 = [cos(phi),sin(phi),0;
%       -sin(phi),cos(phi),0;
%       0,0,1];
% %assume 3-2-1 sequence
% R_321 = R_1*R_2*R_3;
%%from the definition of quaternion matrix form and the 321DCM we can solve for angles
%quaternion vextor defined as= i+ j+ k +w
%%
    a = 2*(q(2)*q(3)+q(1)*q(4));
    b = -q(1)^2-q(2)^2+q(3)^2+q(4)^2;
    YAW = atan2(a,b);%psi
    YAW = radtodeg(YAW);
    
    cc = -2*(q(1)*q(3)-q(2)*q(4));
    %PITCH = asin(cc);%theta
    PITCH = atan2((-2*(q(1)*q(3)-q(2)*q(4))),sqrt(1-(-2*(q(1)*q(3)-q(2)*q(4)))^2));
    PITCH = radtodeg(PITCH);
    
    d = 2*(q(1)*q(2)+q(3)*q(4));
    e = q(1)^2-q(2)^2-q(3)^2+q(4)^2;
    ROLL = atan2(d,e);%phi
    ROLL = radtodeg(ROLL);
    
  zyx=[YAW,PITCH,ROLL]'; %psi,theta,phi in degrees
end

