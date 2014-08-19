%==============================================================================
% (c) Jan Modersitzki 2010/12/28, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% Tutorial for FAIR: 2D interpolation and transformations
%
% - load data                  (setupUSData)
% - generate holomorphic map
% - interpolate                (linearInter) 
% - and visualize              (viewImage2D)
%==============================================================================

setupUSData; 
box = (omega(2:2:end)-omega(1:2:end));
xc = reshape(getCellCenteredGrid(omega,m),[],2);
yc = [(box(1)*((1 - 0.9*xc(:,2)/box(2)).*cos(pi*(1-xc(:,1)/box(1)))/2 + 0.5))
  (box(2)*(1-(1 - 0.9*xc(:,2)/box(2)).*sin(pi*(1-xc(:,1)/box(1)))))];

Tc = linearInter(dataT,omega,yc);
figure(2); viewImage2D(Tc,omega,m,'colormap','gray(256)'); 