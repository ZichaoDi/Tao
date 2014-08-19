% (c) Jan Modersitzki 2010/12/23, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% FAIR-BIG-Tutorials: INTERPOLATION
%
% see also
%   E3_1D_basics, 1D basic interpolation example
%   E3_1D_scale            - 1D multi-scale example
%   E3_1D_derivatives      - 1D check the derivative example
%%                       
%   E3_Hands_ij2xy         - 2D data example, ij <-> xy, multi-level
%   E3_viewImage           - 2D visualize data
%   E3_setupHandData       - load data (with landmarks) and visualize
%%
%   E3_2D_basics           - 2D basic interpolation example
%   E3_2D_scale            - 2D multi-scale example
%   E3_2D_generic          - 2D high-res and low-res representation
%   E3_2D_US_trafo         - 2D US, rotate image
%   E3_2D_derivative       - 2D check the derivative example

clc, clear, close all, help(mfilename);

jobs = {
%   'E3_1D_basics','1D basic interpolation example',...
%   'E3_1D_scale','1D multi-scale example',...
%   'E3_1D_derivatives','1D check the derivative example',...
%   'E3_2D_basics','2D basic interpolation example',...
%   'E3_Hands_ij2xy','2D data example, ij <-> xy, multi-level',...
%   'E2_viewImage','2D visualize data',...
%   'E2_setupHandData','load data (with landmarks) and visualize',...
%   'E3_2D_scale','2D multi-scale example',...
%   'E3_2D_generic','2D high-res and low-res representation',...
%   'E4_US_trafo','2D US, rotate image',...
  'E3_2D_derivative','2D check the derivative example'
};

for k=1:length(jobs)/2,
  runThis(jobs{2*k-1},jobs{2*k});
end;

fprintf('\n<%s> done!\n',mfilename); 
