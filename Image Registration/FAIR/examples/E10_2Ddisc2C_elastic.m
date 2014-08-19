%==============================================================================
% (c) Fabian Gigengack, Jan Modersitzki, Lars Ruthotto 
% 2011/04/13, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% Tutorial for FAIR: Initializes images disc "o" and "C"
%
% Runs hyper elastic multi-level image registration with linear elastic
% regularization (matrixFree)
%==============================================================================

clear, close all, help(mfilename), 

setup2Ddisc2CData;

% prepare the plot
FAIRplots('clear')
Dshow = @(T,R,omega,m) viewIP(abs(T-R),omega,m,'colormap',gray(256));
FAIRplots('set','Dshow',Dshow);

% initialize the regularizer for the non-parametric part
alpha       = 1e3;

[reg,regOptn] = regularizer('reset','regularizer','mfElastic',...
  'alpha',alpha);

% finally: run the MultiLevel Non-Parametric Image Registration
[yc,wc,his] = MLIR(MLdata, 'parametric', false, 'minLevel', 3, 'maxLevel', 6,...
    'maxIterNPIR', 20, 'NPIRLS', @Armijo, 'plots',1);