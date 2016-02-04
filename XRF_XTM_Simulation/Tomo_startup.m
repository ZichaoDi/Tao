disp('X-ray simulation startup file...');
more on;
format compact;
warning off;

mcsDir='/homes/wendydi/Documents/Research';
macDir='/Users/Wendydi/Documents/MATLAB';
if (ispc)
  slash = '\';
else
  slash = '/';
end
if(ismac)
addpath([macDir,slash,'Di_MATLABtool']);
else
addpath([mcsDir,slash,'Di_MATLABtool']);
end

PWD = pwd;
addpath_recurse([pwd, slash, 'result']);
path(path,'../TN');
path(path,'../MGOPT');
% path(path,[pwd, slash, 'data']);
addpath_recurse('./data');
addpath_recurse('../SimplerCode');
if(ismac)
addpath_recurse([macDir,slash,'APSdata']);
addpath_recurse([macDir,slash,'Di_MATLABtool']);
else
addpath_recurse([mcsDir,slash,'APSdata']);
addpath_recurse([mcsDir,slash,'Di_MATLABtool']);
addpath_recurse([mcsDir,slash,'optdi/Reconstruction/Writing/multimodal/version1/GenerateFig']);
end
 ADiMat_startup;






