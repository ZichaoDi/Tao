% (c) Jan Modersitzki 2010/12/23, see FAIR.2 and FAIRcopyright.m.
%   <http://www.mic.uni-luebeck.de/people/jan-modersitzki.html>
%
% Plots Results from Landmark Based Registrations
%
% see also E5_linear, E5_quadratic, E5_TPS

show = @(yc,str) viewImage2D(inter(dataT,omega,yc(:)),omega,m,...
  'colormap','bone(256)','title',str);

subplot(1,2,1); cla; % plot T(X) with LM t_j=LM(j,1:2) and y(r_j)=LM(j,5:6)
show(xc,'T(xc)'); hold on;
ph1 = plot(LM(:,1),LM(:,2),'r+',...
  LM(:,5),LM(:,6),'gx',LM(:,[1,5])',LM(:,[2,6])','r-');

subplot(1,2,2); cla; % plot T(Y) with LM r_j=LM(j,1:2)
show(yc,'T(yc)'); hold on;
ph2 = plot(LM(:,3),LM(:,4),'gx','markersize',30','linewidth',3);
set([ph1;ph2],'markersize',10','linewidth',3);