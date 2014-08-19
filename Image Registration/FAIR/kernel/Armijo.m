

function [t,Yt,LSiter,LS] = Armijo(objFctn,Yc,dY,Jc,dJ,varargin)
global H
if nargin == 0, % help and minimal example
    help(mfilename);   GaussNewton;  return;
end;

LSMaxIter   = 10;           % max number of trials
LSreduction = 1e-4;         % slope of line
t = 1;                      % initial step
for k=1:2:length(varargin), % overwrites default parameter
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

L = [];
descent =   dJ * dY;
for LSiter =1:LSMaxIter,
  Yt = Yc + t*dY;       % compute test value Yt
  Jt = objFctn(Yt);      % evalute objective function
 LS = (Jt<Jc + t*LSreduction*descent); % compare
%   LS = (Jt<Jc);
  if LS, break; end;    % success, return
  L = [L;t,Jt];
  t = t/2;          % reduce t
  
end;
if(~LS)
fprintf('descent direction = %e \n',descent);
vec=10.^-(16:-0.5:1);
Ftest=[];
Ftest1=[];
Ftest2=[];
for i=1:length(vec)
    Ft=Yc + vec(i)*dY;       % compute test value Yt
  Ftest(i) = objFctn(Ft);
  Ftest1(i)=Jc+vec(i)*descent;
  Ftest2(i)=Jc+vec(i)*dY'*H*dY;
end
figure(112)
semilogx(vec,Ftest,'b*--',vec,Ftest1,'r.--',vec,Ftest2,'go--')
title('blue: f(x+hv); red: f(x)+h*df*v; green: f(x)+h*v*H*v;');
end
if LS, return; end;      % we are fine
fprintf('Line Search failed - break\n');
t = 0; Yt = Yc;        % take no action

