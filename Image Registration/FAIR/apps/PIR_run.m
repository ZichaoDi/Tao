clc, clear, close all, help(mfilename);
rand('state',0);
jobs = {
 %   'E2_SSDvsRigid', 'SSD versus rigid (Lena)',...
   'E1_PIR_GN','PIR, Lena, rigid, SSD, Gauss-Newton',
};

for k=1:length(jobs)/2,
  runThis(jobs{2*k-1},jobs{2*k});
end;

fprintf('\n<%s> done!\n',mfilename); 
 
