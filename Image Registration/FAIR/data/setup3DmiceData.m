% ==================================================================================
% (c) Fabian Gigengack and Lars Ruthotto 2011/02/08, see FAIR.2 and FAIRcopyright.m.
% (http://www.uni-muenster.de/EIMI/ ,  http://www.mic.uni-luebeck.de/)
% 
% Cardiac gated [F-18]FDG PET data of a mouse heart acquired with the
% quadHIDAC small animal PET at the European Institute for Molecular
% Imaging (EIMI), University of Muenster, Germany.
% 
% see also data/contents.m

example = '3D-mice';
checkSetupDataFile; if OK, return; end;

% set view options and interpolation options
[viewer,viewOptn] = viewImage('reset','viewImage','imgmontage',...
	'colormap','gray(256)','direction','-zyx');

% setup interpolation scheme
intOptn  = {'inter','splineInterMex','regularizer','moments','theta',.01};

% setup  transformation used in the parametric part
traOptn  = {'trafo','affine3D'};

% setup distance measure
disOptn  = {'distance','SSD'};

% % initialize the regularizer for the non-parametric part
regOptn = {'regularizer','mbHyperElastic','alpha',1,'alphaLength',1,'alphaArea',.1,'alphaVolume',2};

FAIRmessage(mfilename)
load('mice3D');
dataT = double(dataT);
dataR = double(dataR);
MLdata = getMultilevel({dataT,dataR},omega,m,'fig',2);
save(outfile,'dataT','dataR','omega','m','MLdata');
save(outfile,'-append','viewOptn','intOptn','traOptn','disOptn','regOptn');
checkSetupDataFile;

xc       = getCellCenteredGrid(omega,m);
viewData = @(I) viewImage(inter(I,omega,xc),omega,m);

FAIRfigure(1,'figname',mfilename); clf;
subplot(1,2,1); viewData(dataT); title('template');
subplot(1,2,2); viewData(dataR); title('reference');