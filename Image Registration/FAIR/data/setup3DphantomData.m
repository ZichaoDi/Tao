% ==================================================================================
% (c) Lars Ruthotto 2011/02/08, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/
%
% ==================================================================================

example = '3Dphantom';
checkSetupDataFile; if OK, return; end

FAIRmessage(mfilename)
load('phantom3D')

dataT = double(dataT);
dataR = double(dataR);
% setup image viewer,
[viewer,viewOptn] = viewImage('reset','viewImage','imgmontage',...
  'colormap','gray(256)','direction','-zyx');

FAIRfigure(1); clf;
subplot(1,2,1); viewImage(dataT,omega,m,'title','T1');
subplot(1,2,2); viewImage(dataR,omega,m,'title','T2');

MLdata = getMultilevel({dataT,dataR},omega,m,'fig',2,'restrictdim',[1 1 1]);

% setup interpolation scheme
intOptn   = {'inter','splineInterMex','regularizer','moments','theta',1e-1};

% setup  transformation used in the parametric part
traOptn = {'trafo','affine3Dsparse'};

% setup distance measure
disOptn = {'distance','SSD'};

% initialize the regularizer for the non-parametric part
regOptn = {'regularizer','mbDiffusionEPI','alpha',1,'beta',0.1,...
    'distortionDirections',[0,0;1,-1;0,0]};

% save to outfile
save(outfile,'dataT','dataR','omega','m',...
  'MLdata','viewOptn','intOptn','traOptn','disOptn','regOptn');
checkSetupDataFile

