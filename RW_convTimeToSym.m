% Used by "CalculateAandBMatrix.m" to change the format of the
% strings to do substitutions


function ret = RW_convTimeToSym(func)

ret = char(func);

% x y z dx dy dz d2x d2y d2z
% phi theta psu dphi dtheta dpsu

ret = strrep(ret, 'diff(x(t), t, t)', 'd2x');
ret = strrep(ret, 'diff(y(t), t, t)', 'd2y');
ret = strrep(ret, 'diff(z(t), t, t)', 'd2z');

ret = strrep(ret, 'diff(x(t), t)', 'dx');
ret = strrep(ret, 'diff(y(t), t)', 'dy');
ret = strrep(ret, 'diff(z(t), t)', 'dz');  

ret = strrep(ret, 'x(t)', 'x');
ret = strrep(ret, 'y(t)', 'y');
ret = strrep(ret, 'z(t)', 'z');  


ret = strrep(ret, 'diff(wx, t)', 'dwx');
ret = strrep(ret, 'diff(wy, t)', 'dwy');
ret = strrep(ret, 'diff(wz, t)', 'dwz'); 

ret = strrep(ret, 'wx(t)', 'wx');
ret = strrep(ret, 'wy(t)', 'wy');
ret = strrep(ret, 'wz(t)', 'wz');  

ret = strrep(ret, 'diff(q1(t), t)', 'dq1');
ret = strrep(ret, 'diff(q2(t), t)', 'dq2');
ret = strrep(ret, 'diff(q3(t), t)', 'dq3'); 
ret = strrep(ret, 'diff(q4(t), t)', 'dq4');  

ret = strrep(ret, 'q1(t)', 'q1');
ret = strrep(ret, 'q2(t)', 'q2');
ret = strrep(ret, 'q3(t)', 'q3');
ret = strrep(ret, 'q4(t)', 'q4');

ret = strrep(ret, 'diff(wa(t), t)', 'dwa');
ret = strrep(ret, 'diff(wb(t), t)', 'dwb');
ret = strrep(ret, 'diff(wc(t), t)', 'dwc'); 

ret = strrep(ret, 'wa(t)', 'wa');
ret = strrep(ret, 'wb(t)', 'wb');
ret = strrep(ret, 'wc(t)', 'wc');  

ret = strrep(ret, ', t == 0..1/100), t == 0..1/100)','');
ret = strrep(ret, ', t == 0..1/100)','');    
ret = strrep(ret, 'int(int(',''); 
ret = strrep(ret, 'int(',''); 


