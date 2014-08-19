example = 'CircleToC';
checkSetupDataFile; if OK, return; end

% set view options and interpolation options
[viewer,viewOptn] = viewImage('reset','viewImage','viewImage2D',...
  'colormap',flipud(gray(256)));

% setup interpolation scheme
intOptn = {'inter','splineInter','regularizer','moments','theta',1};

% setup  transformation used in the parametric part
traOptn = {'trafo','affine2D'};

% setup distance measure
disOptn = {'distance','SSD'};

% initialize the regularizer for the non-parametric part
[reg,regOptn] = regularizer('reset','regularizer', ...
    'mbHyperElastic','alpha',1,'alphaLength',1000,'alphaArea',0,'alphaVolume',500);

FAIRmessage(mfilename)
get2Ddata(outfile,'circle.jpg','c.jpg', ...
    'omega',[0,1,0,1],'m',[128,128]);
save(outfile,'-append','viewOptn','intOptn','traOptn','disOptn','regOptn');

checkSetupDataFile;