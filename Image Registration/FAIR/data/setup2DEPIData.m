% =============================================================================
% (c) Lars Ruthotto 2011/02/08, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/
%
% loads or generates 2D slices of spin-echo echo planar MR images. Both
% images are affected by geometrical deformations and intensity modulations
% in opposite directions.
%
% Data supplied by Harald Kugel, Institute of Clinical Radiology, Univesity
% Hospital Muenster, Germany.
%
% Data originates from
%  @inproceedings{2010-SPIE-ORKSFW,
%	Author = {Olesch, Janine and Ruthotto, Lars and Kugel, Harald and Skare, Stefan and Fischer, Bernd and Wolters, Carsten H.},
%	Editor = {Benoit M. Dawant and David R. Haynor},
%	Journal = {Medical Imaging 2010: Image Processing},
%	Location = {San Diego, California, USA},
%	Number = {1},
%	Numpages = {8},
%	Pages = {76230K},
%	Publisher = {SPIE},
%	Title = {A Variational Approach for the Correction of Field-inhomogeneities in EPI Sequences},
%	Volume = {7623},
%	Year = {2010}
%  }
% 
% see also apps/EPI/contents
%==============================================================================

example = 'EPIslice';
checkSetupDataFile; if OK , return;  end;

% set view options and interpolation options
[viewer,viewOptn] = viewImage('reset','viewImage','viewImage2Dsc','colormap','gray(256)');

% setup interpolation scheme
intOptn  = {'inter','splineInterMex','regularizer','moments','theta',1e-1};

% setup  transformation used in the parametric part
traOptn  = {'trafo','affine2D'};

% setup distance measure
disOptn  = {'distance','SSD'};

% % initialize the regularizer for the non-parametric part
regOptn = {'regularizer','mbDiffusionEPI','alpha',1,'beta',.1,...
    'distortionDirections',[0,0;1,-1]};

FAIRmessage(mfilename)
get2Ddata(outfile,'EPIslice-T.jpg','EPIslice-R.jpg',...
  'omega',[0,239.7126,0,239.0557],'m',[256,256]);
save(outfile,'-append','viewOptn','intOptn','traOptn','disOptn','regOptn');

checkSetupDataFile;