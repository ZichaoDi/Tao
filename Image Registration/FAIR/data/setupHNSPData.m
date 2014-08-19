% (c) Jan Modersitzki 2010/12/25, see FAIR.2 and FAIRcopyright.m.
% Flexible Algorithms for Image Registration, SIAM, 2009
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
% 
% Images from a histologial serial sectioning, 
% data courtesy Oliver Schmitt, Institute of Anatomy, University of Rostock, Germany
%
% see also get2Ddata and data/contents.m

example = 'HNSP';
checkSetupDataFile; if OK, return; end;

% set view options and interpolation options
viewOptn = {'viewImage','viewImage2D','colormap','gray(256)'};
viewImage('reset',viewOptn{:});
FAIRmessage(mfilename)
get2Ddata(outfile,'HNSP-T.jpg','HNSP-R.jpg','omega',[0,2,0,1],'m',[512,256]);
load(outfile);