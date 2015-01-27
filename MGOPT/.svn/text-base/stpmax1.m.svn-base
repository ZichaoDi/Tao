function alpha = stpmax1 (v, e, v_low, v_up)
global current_n
%------------------------------------------------
% compute the maximum allowable step length
%------------------------------------------------
% Usage: alpha = stpmax1 (v, e, v_low, v_up)
%------------------------------------------------
al = 1e20; % initialize in case "if" statements
au = 1e20; % are not executed
alpha=zeros(size(v));
indl = find(e < 0);
if (length(indl) > 0);
   tl   = v_low(indl) - v(indl);
   %alpha(indl)=tl./e(indl);
   al   = min(tl./e(indl));

end;
%------------------------------------------------
indu = find(e > 0 );
if (length(indu) > 0);
   tu   = v_up (indu) - v(indu);
   %alpha(indu)=tu./e(indu);
   au   = min(tu./e(indu));

end;
%------------------------------------------------
alpha  = min([1 al au]);
alpha=1;

% figure(444);plot(1:current_n,e+v,'r.-');title('new v with step length=1');
% figure(333);
% subplot(4,1,1);plot(1:current_n,e,'r.-');title('search-direction');
% subplot(4,1,2);plot(1:current_n,v,'r.-');title('current_v');
%    nn=length(tl);
%    mm=length(tu);
%    subplot(4,1,4);plot(1:mm,tu,'r.-',1:mm,tu./e(indu),'b.-');title('tu(red);au(blue)');
%    subplot(4,1,3);plot(1:nn,tl,'r.-',1:nn,tl./e(indl),'b.-'); title('tl(red);al(blue)');