% =========================================================================
% (c) Jan Modersitzki 2010/12/27, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% Tutorial for FAIR: REGULARIZER
%
% E8_forces:           elastic responce to a force field
% E8_matrices:         creating B, matrix based an matrix free version
% E8_setup:            regularizarion setup
% E8_checkOperations:  regularizer: multigrid example
% =========================================================================

clc, clear, close all, help(mfilename);

% runThis('E8_forces',   'elastic response to a force field');
 runThis('E8_forces_curvature',   'curvature response to a force field');
% runThis('E8_matrices', 'creating elastic B, matrix based an matrix free version');
% runThis('E8_matrices_curvature', 'creating curvature B, matrix based an matrix free version');
% runThis('E8_setup',    'elastic regularizer setup');
% runThis('E8_setup_curvature',    'curvature regularizer setup');
% runThis('E8_checkOperations',    'regularizer: multigrid example');

fprintf('\n<%s> done!\n',mfilename); 