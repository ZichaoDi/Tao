% (c) Jan Modersitzki 2010/12/25, see FAIR.2 and FAIRcopyright.m.
% Flexible Algorithms for Image Registration, SIAM, 2009
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
% 
% Initializing MRI slices of a human head (T1/T2).
%
% see also get2Ddata and data/contents.m

example = 'MRI-head';
checkSetupDataFile; if OK, return; end;

% set view options and interpolation options
viewOptn = {'viewImage','viewImage2D','colormap','bone(256)'};
viewImage('reset',viewOptn{:});
FAIRmessage(mfilename)
get2Ddata(outfile,'MRIhead-T.jpg','MRIhead-R.jpg','omegaT',[0,20,0,20],'omegaR',[0,20,0,20]);
load(outfile);
