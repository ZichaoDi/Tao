%==============================================================================
% (c) Jan Modersitzki 2010/12/28, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% Tutorial for FAIR: translation of an US image
%
% - load data                  (setupUSData)
% - transform                  (plain translation)
% - interpolate                (linearInter) 
% - and visualize              (viewImage2D)
%==============================================================================

setupUSData; 
wc = [-50;0]; 
xc = reshape(getCellCenteredGrid(omega,m),[],2);  
yc = [(xc(:,1) + wc(1));(xc(:,2) + wc(2))];
Tc = linearInter(dataT,omega,yc);
figure(2); clf; viewImage2D(Tc,omega,m,'colormap','gray(256)'); 
  