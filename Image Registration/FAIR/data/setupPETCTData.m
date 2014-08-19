% (c) Jan Modersitzki 2010/12/25, see FAIR.2 and FAIRcopyright.m.
% Flexible Algorithms for Image Registration, SIAM, 2009
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
% 
% PET-CT data of a human thorax. Original data from 
% 	@article{ShekharEtAl2005,
% 	 author = {Raj Shekhar and Vivek Walimbe and Shanker Raja and Vladimir Zagrodsky
% 	           and Mangesh Kanvinde and Guiyun Wu and Bohdan Bybel},
% 	 title = {Automated 3-Dimensional Elastic Registration of Whole-Body {PET}
% 	          and {CT} from Separate or Combined Scanners},
% 	 journal = {J. of Nuclear Medicine},
% 	 volume = {46},
% 	 number = {9},
% 	 year = {2005},
% 	 pages = {1488--1496},
% 	}
%
% see also get2Ddata and data/contents.m

example = 'PETCT';
checkSetupDataFile; if OK, return; end;

% set view options and interpolation options
viewOptn = {'viewImage','viewImage2D','colormap','gray(256)'};
viewImage('reset',viewOptn{:});
FAIRmessage(mfilename)
get2Ddata(outfile,'PET-CT-PET.jpg','PET-CT-CT.jpg',...
  'omegaT',[0,50,0,50],'omegaR',[0,50,0,50],'m',[128,128]);
load(outfile);

