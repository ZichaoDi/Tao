% (c) Jan Modersitzki 2010/12/25, see FAIR.2 and FAIRcopyright.m.
% Flexible Algorithms for Image Registration, SIAM, 2009
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
% 
% MRI data of a human brain, 
% data courtesy Ron Kikinis, Brigham & Women's Hospital, Boston, USA
%
% Data originates from Ron Kikinis, 

example = '3Dbrain';
checkSetupDataFile; if OK,   return; end;

FAIRmessage(mfilename)

load('brain3D')

% setup image viewer,
[viewer,viewOptn] = viewImage('reset','viewImage','imgmontage',...
  'colormap','gray(256)','direction','-zyx');

FAIRfigure(1); clf;
subplot(1,2,1); viewImage(dataT,omega,m);
subplot(1,2,2); viewImage(dataR,omega,m);

MLdata = getMultilevel({dataT,dataR},omega,m,'fig',2);

% setup interpolation scheme
intOptn   = {'inter','splineInterMex'};

% setup  transformation used in the parametric part
traOptn = {'trafo','affine3Dsparse'};

% setup distance measure
disOptn = {'distance','SSD'};

% initialize the regularizer for the non-parametric part
regOptn = {'regularizer','mfElastic','alpha',1e3,'mu',1,'lambda',0};

% save to outfile
save(outfile,'dataT','dataR','omega','m',...
  'MLdata','viewOptn','intOptn','traOptn','disOptn','regOptn');
checkSetupDataFile