disp('X-ray simulation startup file...');
more on;
format compact;
warning off;

if (ispc)
  slash = '\';
else
  slash = '/';
end

PWD = pwd;
path(path,[pwd, slash, 'data']);
path(path,[pwd, slash, 'result']);
path(path,'../TN');





