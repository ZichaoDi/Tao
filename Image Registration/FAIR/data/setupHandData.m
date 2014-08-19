% (c) Jan Modersitzki 2010/12/25, see FAIR.2 and FAIRcopyright.m.
% Flexible Algorithms for Image Registration, SIAM, 2009
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
% 
% Initializing xray data of human hands and adds manually chosen landmarks.
% Data originates from
%  @article{Amit1994,
%    author = {Yali Amit},
%     title = {A nonlinear variational problem for image matching},
%      year = {1994},
%   journal = {SIAM J. Sci. Comput.},
%    volume = {15},
%    number = {1},
%     pages = {207--224},
%  }
% 
% see also get2Ddata and data/contents.m

example = 'Hands';
checkSetupDataFile; if OK, return; end;

% set view options and interpolation options
viewOptn = {'viewImage','viewImage2D','colormap','bone(256)'};
viewImage('reset',viewOptn{:});

FAIRmessage(mfilename)

get2Ddata(outfile,'hands-T.jpg','hands-R.jpg',...
  'omega',[0,20,0,25],'m',[128,128]);
% get2Ddata(outfile,'LenaCropped.tiff','LenaReference.tiff',...
%   'omega',[0,50,0,50],'m',[512,512]);

% for this data landmarks have been identified
LM    = [
  5.5841   17.2664    2.6807   12.7797
  10.7243   21.6121    7.2028   19.6795
  13.2477   21.6121   10.1865   20.8916
  15.2570   19.2290   12.5175   20.0991
  15.8645   15.1636   14.3357   16.7424
  5.3972    8.1075    7.9953    6.3462
  7.5000    5.9579   11.8648    5.6469
  ];

% save to outfile
save(outfile,'-append','LM','example');
load(outfile);