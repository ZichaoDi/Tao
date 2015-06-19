disp('X-ray simulation startup file...');
more on;
format compact;
warning off;

if (ispc)
  slash = '\';
else
  slash = '/';
end
if(ismac)
addpath('/Users/Wendydi/Documents/MATLAB/Di_MATLABtool');
else
addpath('/homes/wendydi/Documents/Research/Di_MATLABtool');
end

PWD = pwd;
path(path,[pwd, slash, 'result']);
path(path,'../TN');
path(path,'../MGOPT');
% path(path,[pwd, slash, 'data']);
addpath_recurse('./data');
addpath_recurse('./result');
addpath_recurse('../SimplerCode');
if(ismac)
addpath_recurse('/Users/Wendydi/Documents/MATLAB/APSdata');
addpath_recurse('/Users/Wendydi/Documents/MATLAB/Di_MATLABtool');
else
addpath_recurse('/homes/wendydi/Documents/Research/Di_MATLABtool');
addpath_recurse('/homes/wendydi/Documents/Research/APSdata');
end
 ADiMat_startup;






