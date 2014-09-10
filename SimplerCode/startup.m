disp('Reconstruction startup file...');
more on;
format compact;
warning off;

%path(path,'C:\Documents and Settings\buckaroo\Desktop\mgopt_04_06\MGOPT')
%path(path,'C:\Documents and Settings\buckaroo\Desktop\mgopt_04_06\TN')
%path(path,'C:\Documents and Settings\buckaroo\Desktop\mgopt_04_06\NSOLA')

if (ispc)
  slash = '\';
else
  slash = '/';
end

PWD = pwd;
path(path,[pwd, slash, 'imagine']);
path(path,[pwd, slash, 'data']);
path(path,[pwd, slash, 'data',slash,'attData']);





