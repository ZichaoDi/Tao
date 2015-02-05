function [ipivot1, flast1] = modz (x, p, ipivot, low, up, flast, f, alpha)
%---------------------------------------------------------------------
% update the constraint matrix if a new constraint is encountered
%---------------------------------------------------------------------
indl = [1:length(p)]';%find(ipivot == 0 & p < 0);%
if (length(indl) > 0);
   toll = 10 * eps^(1/2) * (abs(low(indl)) + ones(size(indl)));
   hitl = find(x(indl)-low(indl) <= toll);
   if (length(hitl)>0);
      flast = f;
      ipivot(indl(hitl)) = -1;
   end;
   indDrop=find(ipivot==-1 & x-low > 10 * eps^(1/2));
   if(~isempty(indDrop))
       ipivot(indDrop)=0;
   end
end;
%---------------------------------------------------------------------
indu = find(ipivot == 0 & p > 0);%[1:length(p)]';%
if (length(indu) > 0);
   tolu = 10 * eps * (abs( up(indu)) + ones(size(indu)));
   hitu = find(up(indu)-x(indu)  <= tolu);
   if (length(hitu)>0);
      flast = f;
      ipivot(indu(hitu)) = 1;
   end;
end;
%---------------------------------------------------------------------
flast1  = flast;
ipivot1 = ipivot;
