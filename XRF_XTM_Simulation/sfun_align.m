function [F,G,f]=sfun_align(x,numThetan,thetan,w)
f=zeros(2*numThetan,1);
f(1:numThetan)=x(3)+x(1)*sin(thetan'*pi/180+x(2))-w(:,1);
f(numThetan+1:2*numThetan)=x(6)+x(4)*sin(thetan'*pi/180+x(5))-w(:,2);
G=zeros(6,1);
G(1)=2*sin(thetan.*pi/180+x(2))*f(1:numThetan);
G(2)=2*x(1)*cos(thetan.*pi/180+x(2))*f(1:numThetan);
G(3)=2*sum(f(1:numThetan));
G(4)=2*sin(thetan.*pi/180+x(5))*f(numThetan+1:end);
G(5)=2*x(4)*cos(thetan.*pi/180+x(5))*f(1+numThetan:end);
G(6)=2*sum(f(1+numThetan:end));

F=sum(f.^2);
%%============================= Two Wires
% F=zeros(numThetan*2,1);
% F(1:numThetan)=x(5:end)+x(1)*sin(thetan'*pi/180+x(2))-w(:,1);
% F(numThetan+1:end)=x(5:end)+x(3)*sin(thetan'*pi/180+x(4))-w(:,2);
% G=zeros(numThetan+4,1);
% G(1)=2*sin(thetan*pi/180+x(2))*F(1:numThetan);
% G(2)=2*x(1)*cos(thetan*pi/180+x(2))*F(1:numThetan);
% G(3)=2*sin(thetan*pi/180+x(4))*F(numThetan+1:end);
% G(4)=2*x(3)*cos(thetan*pi/180+x(3))*F(numThetan+1:end);
% G(5:end)=2*F(1:numThetan)+2*F(numThetan+1:end);
% F=sum(F.^2);

% F=zeros(numThetan*2,1);
% for n=1:numThetan,
%     F(n)=x(1)*sin(thetan(n)*pi/180+x(2))-w(n,1);
%     F(numThetan+n)=x(3)*sin(thetan(n)*pi/180+x(4))-w(n,2);
% end
% G=zeros(4,1);
% G(1)=2*sin(thetan*pi/180+x(2))*F(1:numThetan);
% G(2)=2*x(1)*cos(thetan*pi/180+x(2))*F(1:numThetan);
% G(3)=2*sin(thetan*pi/180+x(4))*F(numThetan+1:end);
% G(4)=2*x(3)*cos(thetan*pi/180+x(3))*F(numThetan+1:end);
% F=sum(F.^2);
