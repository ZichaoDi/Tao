% ==================================================================================
% (c) Jan Modersitzki and Fabian Gigengack and Lars Ruthotto 2010/12/25
% see FAIR.2 and FAIRcopyright.m.
% Flexible Algorithms for Image Registration, SIAM, 2009
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
% http://www.uni-muenster.de/EIMI/
% http://www.mic.uni-luebeck.de/

example = 'MPsquare';
checkSetupDataFile; if OK&&0, return; end

% set view options and interpolation options
[viewer,viewOptn] = viewImage('reset','viewImage','viewImage2D','colormap','gray(256)');

% setup interpolation scheme
intOptn = {'inter','splineInterMex','regularizer','moments','theta',1e-2};

% setup  transformation used in the parametric part
traOptn = {'trafo','affine2D'};

% setup distance measure
disOptn = {'distance','SSD'};

% initialize the regularizer for the non-parametric part
[reg,regOptn] = regularizer('reset','regularizer', ...
    'mfHyperElastic','alpha',1,'alphaLength',100,'alphaArea',0,'alphaVolume',100);

FAIRmessage(mfilename)
get2Ddata(outfile,'MPSquare-T.jpg','MPSquare-R.jpg', ...
    'omega',[0 64 0 64],'m',[64 64]);
save(outfile,'-append','viewOptn','intOptn','traOptn','disOptn','regOptn');

checkSetupDataFile;