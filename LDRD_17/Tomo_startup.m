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
addpath_recurse(['./result']);
path(path,'../TN');
path(path,'../MGOPT');
addpath_recurse('../XRF_XTM_Simulation/data');
addpath_recurse('../Recon_onfly/result');
addpath_recurse('../SimplerCode');
if(ismac)
addpath_recurse([macDir,slash,'APSdata']);
addpath_recurse([macDir,slash,'Di_MATLABtool']);
else
addpath_recurse([mcsDir,slash,'APS']);
addpath_recurse([mcsDir,slash,'Di_MATLABtool']);
addpath_recurse([mcsDir,slash,'multigrid']);
addpath_recurse([mcsDir,slash,'optdi/Reconstruction/Writing/multimodal/version1/GenerateFig']);
% addpath([mcsDir,slash,'Di_MATLABtool/dip/common/dipimage']);
% dip_initialise;
% dipsetpref('imagefilepath' , '/images');
end
% ADiMat_startup;






