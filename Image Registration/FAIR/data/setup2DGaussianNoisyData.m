% ==================================================================================
% (c) Fabian Gigengack and Lars Ruthotto 2011/02/08, see FAIR.2 and FAIRcopyright.m.
% (http://www.uni-muenster.de/EIMI/ ,  http://www.mic.uni-luebeck.de/)

example = 'GaussFunctionsNoisy';
checkSetupDataFile; if OK, return; end

% set view options and interpolation options
[viewer,viewOptn] = viewImage('reset','viewImage','viewImage2D','colormap','gray(256)');

% setup interpolation scheme
intOptn = {'inter','splineInter','regularizer','moments','theta',10};

% setup  transformation used in the parametric part
traOptn = {'trafo','affine2D'};

% setup distance measure
disOptn = {'distance','SSD'};

% initialize the regularizer for the non-parametric part
[reg,regOptn] = regularizer('reset','regularizer', ...
    'mfHyperElastic','alpha',1,'alphaLength',1000,'alphaArea',0,'alphaVolume',1000);

FAIRmessage(mfilename)
get2Ddata(outfile,'Gauss-T-noisy.jpg','Gauss-R-noisy.jpg', ...
    'omega',[-5,5,-5,5],'m',[256,256]);
save(outfile,'-append','viewOptn','intOptn','traOptn','disOptn','regOptn');

checkSetupDataFile;