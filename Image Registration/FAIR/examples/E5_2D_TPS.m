%==============================================================================
% (c) Jan Modersitzki 2010/12/28, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% Tutorial for FAIR:  Landmark Based Registration, TPS transformation (various theta's)
%
% - load data (see setupHandData)
% - setup  viewer (viewImage2D), interpolator (splineInter), 
% - setup landmarks (LM)
% - run TPS, various theta
%==============================================================================

clear, close all, help(mfilename)

%% setup hand data
setupHandData

if runThis(mfilename,'set new landmarks ? ','test',0),
  [LM,fig] = getLandmarks(dataT,dataR,omega,m);
  close(fig);
end;

omegaT = omega(1,:);
omegaR = omega(end,:);
xT = getCellCenteredGrid(omegaT,m);
xR = getCellCenteredGrid(omegaR,m);
Tc = inter(dataT,omegaT,xT);
Rc = inter(dataR,omegaR,xR);

%% visualize data
FAIRfigure(1,'figname',mfilename); clf; 
subplot(1,3,1); viewImage(Tc,omegaT,m); hold on;
ph = plotLM(LM(:,1:2),'numbering','on','color','r');
set(ph,'linewidth',2,'markersize',20);
title(sprintf('%s','T&LM'),'fontsize',20);

subplot(1,3,2); viewImage(Rc,omegaR,m); hold on;
ph = plotLM(LM(:,3:4),'numbering','on','color','g','marker','+');
set(ph,'linewidth',2,'markersize',20);
title(sprintf('%s','R&LM'),'fontsize',20);

%% run over different theta's
mh = [];
theta = [1e3,1e2,1e1,0];
for j=1:length(theta),
  [yc,LM] = LMreg('TPS',LM(:,1:4),xR,'theta',theta(j));
  TLM = inter(dataT,omegaT,yc);

  set(mh,'visible','off');
  subplot(1,3,3); cla; viewImage(TLM,omegaR,m); hold on;
  ph = plotLM(LM(:,3:4),'numbering','off','color','g','marker','+');
  qh = plotLM(LM(:,7:8),'numbering','off','color','m','marker','x');
  rh = plot(LM(:,[3,7])',LM(:,[4,8])','m-','linewidth',3);
  set([ph;qh;rh],'linewidth',2,'markersize',20);
  title(sprintf('T(Y^{TPS},\\theta=%s)&LM',num2str(theta(j))),'fontsize',20);
  pause;
end;
