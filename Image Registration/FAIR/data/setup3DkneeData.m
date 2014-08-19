% (c) Jan Modersitzki 2010/12/25, see FAIR.2 and FAIRcopyright.m.
% Flexible Algorithms for Image Registration, SIAM, 2009
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
% 
% MRI data of a human knee, image courtesy by 
% Thomas Netsch, Philips Medical Solutions, Hamburg
%
% see also data/contents.m

example = '3Dknee';
checkSetupDataFile; %if OK,   return; end;

FAIRmessage(mfilename)

% setup image viewer,
[viewer,viewOptn] = viewImage('reset','viewImage','imgmontage',...
  'colormap','bone(256)','direction','-zyx');

infile = fullfile(FAIRpath,'data','createPhilipsKnees.mat');
if ~exist(infile,'file'), infile = createPhilipsKnees; end;
load(infile);

dataT = double(dataT);
dataR = double(dataR);

viewK = @(T) volView(T,omega,m,'isovalue',15,'view',[-20,-5],'colormap','bone(256)');

FAIRfigure(1,'color','w','figname',sprintf('%s/template',mfilename)); clf; 
viewK(dataT); hold on; axis off; colormap(gray(128));
FAIRfigure(2,'color','w','figname',sprintf('%s/reference',mfilename)); clf; 
viewK(dataR); hold on; axis off; colormap(gray(128));  

MLdata = getMultilevel({dataT,dataR},omega,m,'fig',3);

% setup interpolation scheme
intOptn   = {'inter','linearInterMex'};

% setup  transformation used in the parametric part
traOptn = {'trafo','affine3Dsparse'};

% setup distance measure
disOptn = {'distance','SSD'};

% initialize the regularizer for the non-parametric part
regOptn = {'regularizer','mfElastic','alpha',500,'mu',1,'lambda',0};

% save to outfile
save(outfile,'dataT','dataR','omega','m',...
  'MLdata','viewOptn','intOptn','traOptn','disOptn','regOptn');
checkSetupDataFile



