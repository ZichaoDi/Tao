disp('X-ray simulation startup file...');
more on;
format compact;
warning off;

if (ispc)
  slash = '\';
else
  slash = '/';
end
addpath('/Users/Wendydi/Documents/MATLAB/Di_MATLABtool');

PWD = pwd;
% path(path,[pwd, slash, 'data']);
path(path,[pwd, slash, 'result']);
path(path,'../TN');
addpath_recurse('./data');
addpath_recurse('../SimplerCode');
addpath_recurse('/Users/Wendydi/Documents/MATLAB/APSdata');
addpath_recurse('/Users/Wendydi/Documents/MATLAB/Di_MATLABtool');
ADiMat_startup;






