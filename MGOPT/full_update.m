function v_h = full_update(v_H)
%--------------------------------------------------
% update from current grid directly to finest grid
%--------------------------------------------------
% Usage: v_h = full_update(v_H)
%--------------------------------------------------
global N current_n
%--------------------------------------------------
n = current_n;
j = find(N==current_n);
%--------------------------------------------------
if (j==1);
   v_h = v_H;
elseif (j==2);
   v_h = update(v_H,0);
else
   for i=j:-1:3;
      v_H = update(v_H,0);
      j = j-1;
      current_n = N(j);
   end;
   v_h = update(v_H,0);
end;
current_n = N(1); 