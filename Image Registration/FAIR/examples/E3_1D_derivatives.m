%==============================================================================
% (c) Jan Modersitzki 2010/12/27, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% Tutorial for FAIR: 1D interpolation, derivatives
%
% - setup 1D data
% - visualize data
% - visualize derivatives
% - compute interpolants Tfctn(x) = spline(TD,omega,x) and check derivative
%==============================================================================

clear, close all, help(mfilename); echo on

%% setup data
fprintf('%s\n','generate data');

%                     x
%                                              T
%             x               x 
% |---x---|---+---|---+---|---+---|---x---|    x
% 0                                       omega

omega = [0,14];         % x data is located in intervall [omega(1),omega(2)]
dataX = 1:2:omega(2);   % data location x(j) and values T(j)
dataT = [0;0;1;4;1;0;0];% the data
fprintf('[dataX;dataT]=\n'); disp(reshape([dataX(:);dataT(:)],[],2)')

% prepare continuous medel
coefT   = getSplineCoefficients(dataT,'dim',1,'out',0);
Tspline = @(x) splineInter(coefT,omega,x);

% output:display
fprintf('%s\n','visualize the data')

% output:plot
figure(1); clf; 
ph = plot(dataX,0*dataT,'k.',dataX,dataT,'b.','markersize',20); hold on; 
title('1D interpolation with derivatives','fontsize',20); 
axis([omega(1)-1,omega(2)+1,min(dataT)-1,max(dataT)+1]);
xFine = getCellCenteredGrid(omega,101);
ph(3) = plot(xFine,Tspline(xFine),'b-','linewidth',3);
lstr = {'location','data','T'}; legend(ph,lstr,1);
pause;

%% play with derivatives: get some interestin points
yc = getCellCenteredGrid(omega,9);
[Tc,dT] = Tspline(yc);

% grep diagonal of dT
dT = full(diag(dT)); dx = 0.35*diff(omega)./length(yc);

for j=1:length(yc),
  qh=plot(yc(j),Tc(j),'m.',yc(j)+dx*[-1,1],Tc(j)+dT(j)*dx*[-1,1],'m-',...
    'markersize',20,'linewidth',2);
end;
lstr{end+1}= 'derivative';
legend([ph;qh(end)],lstr,1)
pause;

% the ultimative derivative check: for MATLAB and FAIR
x = omega(1)+diff(omega)*randn(15,1);
Mfctn = @(x) linearInterMatlab(coefT,omega,x);
checkDerivative(Mfctn,x);
ylabel('MATLAB linear interpolation','fontsize',30);
set(gca,'fontsize',30);

checkDerivative(Tspline,x);
ylabel('FAIR spline interpolation','fontsize',30);
set(gca,'fontsize',30);

echo off