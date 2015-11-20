syms We t
syms Tdx Tdy Tdz
syms Tcx Tcy Tcz
syms Ix Iy Iz
syms q1 q2 q3 q4
%put together vects
wx = sym('wx(t)')
wy = sym('wy(t)')
wz = sym('wz(t)')
W_e_s = [wx; wy; wz]

% Td = [Tdx; Tdy; Tdz]
% Tc = [Tcx; Tcy; Tcz]
% I = diag([Ix, Iy, Iz])
q1 = sym('q1(t)')
q2 = sym('q2(t)')
q3 = sym('q3(t)')
q4 = sym('q4(t)')
q = [q1, q2, q3, q4]
R_s_e = [q(1)^2-q(2)^2-q(3)^2+q(4)^2, 2*(q(1)*q(2)+q(3)*q(4)), 2*(q(1)*q(3)-q(2)*q(4)); 
         2*(q(1)*q(2)-q(3)*q(4)), -q(1)^2+q(2)^2-q(3)^2+q(4)^2, 2*(q(2)*q(3)+q(1)*q(4)); 
         2*(q(1)*q(3)+q(2)*q(4)), 2*(q(2)*q(3)+q(1)*q(4)), -q(1)^2-q(2)^2+q(3)^2+q(4)^4] ;


W_i_e = [0;0;We]
W_i_s = R_s_e*W_i_e+W_e_s;

dW_i_s = diff(W_i_s, t)

s_dW_i_s = convTimeToSym(dW_i_s)

dW = [We*(2*q1*dq3 + 2*q3*dq1 - 2*q2*dq4 - 2*q4*dq2) + dwx; 
      We*(2*q1*dq4 + 2*q2*dq3 + 2*q3*dq2 + 2*q4*dq1) + dwy;
      dwz - We*(2*q1*dq1 + 2*q2*dq2 - 2*q3*dq3 - 4*q4^3*dq4)]