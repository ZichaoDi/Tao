%==============================================================================
% (c) Jan Modersitzki 2011/01/02, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% Tutorial for FAIR: SSD versus rotations, linearInter, HNSP level=8
%
%==============================================================================

clear, close all, help(mfilename);

setupHNSPData; 
inter('set','inter','linearInter'); 
level = 8; omega = MLdata{level}.omega; m = MLdata{level}.m; 
[T,R] = inter('coefficients',MLdata{level}.T,MLdata{level}.R,omega);
xc = getCellCenteredGrid(omega,m); 
Rc = inter(R,omega,xc);

center = (omega(2:2:end)-omega(1:2:end))'/2;
trafo('set','trafo','rotation2D','c',center);

wc = pi/2*linspace(-1,1,101);  dc = zeros(size(wc));
figure(1); clf;
for j=1:length(wc),
  yc = trafo(wc(j),xc);
  Tc = inter(T,omega,yc);
  dc(j) = SSD(Tc,Rc,omega,m);
  viewImage(255-abs(Tc-Rc),omega,m); drawnow; pause(1/60)
end;
figure(2); clf; p1 = plot(wc,dc); 
