% (c) Jan Modersitzki 2010/12/25, see FAIR.2 and FAIRcopyright.m.
% Flexible Algorithms for Image Registration, SIAM, 2009
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
% 
% Academical 3D images: a box and its translate.
%
% see also data/contents.m


example = '3Dbox';
checkSetupDataFile; 
if OK, 
  viewImage('set','viewImage','imgmontage','direction','-zyx','colormap','bone(256)');
  return; 
end;

FAIRmessage(mfilename)

[viewer,viewOptn] = viewImage('reset','viewImage','imgmontage','colormap','bone(256)');

omega = [0,1,0,1,0,1];
m     = [64,64,64];
xc    = getCellCenteredGrid(omega,m);
xc    = reshape(xc,[],3);
dataT = (abs(xc(:,1)-0.5) < 0.25) & (abs(xc(:,2)-0.5) < 0.2)& (abs(xc(:,3)-0.5) < 0.22);
xc    = xc(:);
dataT = 200*reshape(dataT,m);
wc    = reshape([eye(3),[-0.1;0;0]]',[],1);
yc    = affine3D(wc,xc);
dataR = reshape(linearInter(dataT,omega,yc),m);

FAIRfigure(1); clf;
subplot(1,2,1); viewImage(dataT,omega,m);
subplot(1,2,2); viewImage(dataR,omega,m);

MLdata = getMultilevel({dataT,dataR},omega,m,'fig',2);

% setup interpolation scheme
intOptn   = {'inter','linearInter'};

% setup  transformation used in the parametric part
traOptn = {'trafo','affine3D'};

% setup distance measure
disOptn = {'distance','SSD'};

% initialize the regularizer for the non-parametric part
regOptn = {'regularizer','mfElastic','alpha',500,'mu',1,'lambda',0};

% save to outfile
save(outfile,'dataT','dataR','omega','m',...
  'MLdata','viewOptn','intOptn','traOptn','disOptn','regOptn');
checkSetupDataFile