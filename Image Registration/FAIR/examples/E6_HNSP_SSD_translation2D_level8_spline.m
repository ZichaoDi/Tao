%==============================================================================
% (c) Jan Modersitzki 2011/01/02, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% Tutorial for FAIR: SSD versus translations, splineInter, HNSP level=8
%
%==============================================================================

clear, close all, help(mfilename)

setupHNSPData; 
level = 8; m = MLdata{level}.m; 
inter('set','inter','splineInter'); 
[T,R] = inter('coefficients',MLdata{level}.T,MLdata{level}.R,omega);
xc = getCellCenteredGrid(omega,m); Rc = inter(R,omega,xc);

trafo('set','trafo','translation2D');
figure(1); clf;
[w1,w2] = ndgrid(0.2*linspace(-1,1,21),0.2*linspace(-1,1,21));
dc  = zeros(size(w1));
for j=1:numel(dc),
  yc = trafo([w1(j);w2(j)],xc);
  Tc = inter(T,omega,yc);
  dc(j) = SSD(Tc,Rc,omega,m);
  viewImage(Tc,omega,m); pause(1/100)
end;
figure(1); clf; surf(w1,w2,dc); hold on; grid off; contour(w1,w2,dc)
title(sprintf('translation, m=[%d,%d]',m)); view(-135,33);
