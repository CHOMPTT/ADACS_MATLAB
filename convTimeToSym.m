% Convert the time to symbolic variables

function ret = convTimeToSym(func)

ret = char(func);

ret = strrep(ret, 'diff(wx(t), t)', 'dwx');
ret = strrep(ret, 'diff(wy(t), t)', 'dwy');
ret = strrep(ret, 'diff(wz(t), t)', 'dwz'); 

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

ret = strrep(ret, ', t == 0..1/100), t == 0..1/100)','');
ret = strrep(ret, ', t == 0..1/100)','');    
ret = strrep(ret, 'int(int(',''); 
ret = strrep(ret, 'int(',''); 


