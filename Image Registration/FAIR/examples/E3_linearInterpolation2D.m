%==============================================================================
% (c) Jan Modersitzki 2010/12/27, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% Tutorial for FAIR: linear interpolation in 2D
%
%==============================================================================

dataT = flipud([1,2,3,4;1,2,3,4;4,4,4,4])'; 
m     = size(dataT); 
omega = [0,m(1),0,m(2)]; 
omega1 = [0,m(1)+10,0,m(2)+10]; 
M     = {m,10*m};   % two resolutions, coarse and fine
xc = reshape(getCellCenteredGrid(omega,M{1}),[],2);     % coarse resolution
xf = reshape(getCellCenteredGrid(omega,M{2}),[M{2},2]); % fine   resolution
Tc = linearInter(dataT,omega1,xc(:));
Tf = linearInter(dataT,omega1,xf(:));
% clf; ph = plot3(xc(:,1),xc(:,2),Tc(:),'ro'); hold on;
% %qh = surf(xf(:,:,1),xf(:,:,2),reshape(Tf,M{2}));
xc
dataT
image(xc(:,1),xc(:,2),dataT)