% function FAIRmake(task, varargin)
% 
% (c) Fabian Gigengack and Lars Ruthotto 2011/02/04, see FAIR.2 and FAIRcopyright.m.
% http://www.uni-muenster.de/EIMI/
% http://www.mic.uni-luebeck.de/
%
% Compile c-files or delete compiled c-files.
% 
% Input:
% 	task     -               'all' - all files
%                      'distances' - distance files
%                  'interpolation' - interpolation files
%                          'inter'
%                 'regularization' - regularization files
%                    'regularizer'
%                'transformations' - transformations files
%                          'trafo'
%               <single_file_name> - single file
%   varargin - variable input arguments
%
% Variable Arguments:
%    clean   - remove binarys
%    cores   - compile for multi-core processos (openMP)
%    verbose - show MEX infos
% 
% Examples:
%   FAIRmake('all'); % compile all
%   FAIRmake('all','cores',2,'verbose',true); % compile all for 2 cores
%   FAIRmake('inter','clean',true); % delete compiled interpolation files
%   FAIRmake('linearInterMexC.cpp','clean',true); % delete compiled file
%   FAIRmake('linearInterMexC.cpp'); % compile 'linearInterMexC.cpp'

function FAIRmake(task, varargin)

clean   = false; % delete compiled c-files
cores   = 1;     % choose number of cores; if cores>1 openMP is used
verbose = false; % some output

for k=1:2:length(varargin) % overwrites default parameter
    eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end

switch task
    case 'all' % process all c-files
        % distances
        make('rhoSplineC.cpp', clean, cores,verbose);
        make('NGFdotMexC.cpp', clean, cores,verbose);
        make('SSDmexC.cpp', clean, cores,verbose);
        make('NCCmexC.cpp', clean, cores,verbose);
        % interpolation
        make('linearInterMexC.cpp', clean, cores,verbose);
        make('linearInterSmoothMexC.cpp', clean, cores,verbose);
        make('cubicInterMexC.cpp', clean, cores,verbose);
        make('splineInterMexC.cpp', clean, cores,verbose);
        % regularization
        make('geometrymexC.cpp', clean, cores,verbose);
        % transformation
        make('tensorProdC.c', clean, cores,verbose);
    case 'distances' % process distance c-files
        make('rhoSplineC.cpp', clean, cores,verbose);
        make('NGFdotMexC.cpp', clean, cores,verbose);
        make('SSDmexC.cpp', clean, cores,verbose);
        make('NCCmexC.cpp', clean, cores,verbose);
    case {'interpolation', 'inter'} % process interpolation c-files
        % interpolation
        make('linearInterMexC.cpp', clean, cores,verbose);
        make('linearInterSmoothMexC.cpp', clean, cores,verbose);
        make('cubicInterMexC.cpp', clean, cores,verbose);
        make('splineInterMexC.cpp', clean, cores,verbose);
    case {'regularization', 'regularizer'} % process regularizer c-files
        make('geometrymexC.cpp', clean, cores,verbose);
        make('geometrymexCv2.cpp', clean, cores,verbose);
    case {'transformations', 'trafo'} % process transformation c-files
        make('tensorProdC.c', clean, cores,verbose);
    otherwise % process single c-file
        make(task, clean, cores,verbose);
end

function make(filename, clean, cores,verbose)
[pth, name] = fileparts(which(filename));
cpth = pwd; cd(pth) % change to folder with the c-file
[~, maxsize] = computer; % determine if 64bit (maxsize==2^48-1) or 32bit
if maxsize==2^48-1, options = '-largeArrayDims'; else options = ''; end
if verbose>0, options = [ options ' -v ']; end 
if clean % delete existing compiled file
    delete([name '.' mexext])
else % compile file
    if cores>1 % openMP
        setenv('OMP_NUM_THREADS',num2str(cores));
        str = ['mex ' filename ' -O CC=gcc CXX=g++  LD=g++ CFLAGS="\$CFLAGS -fopenmp -ftree-vectorize " CXXFLAGS="\$CXXFLAGS -fopenmp -ftree-vectorize "  LDFLAGS="\$LDFLAGS -fopenmp -ftree-vectorize " ' options];
        if verbose>=0, disp(str); end;
        eval(str);
    else
        setenv('OMP_NUM_THREADS','1');
%        str = ['mex ' filename ' -g  CFLAGS="\$CFLAGS -p -ftree-vectorize " CXXFLAGS="\$CXXFLAGS -p -ftree-vectorize "  LDFLAGS="\$LDFLAGS -p -ftree-vectorize " ' options];
        str = ['mex ' filename '  -O  CFLAGS="\$CFLAGS  -ftree-vectorize " CXXFLAGS="\$CXXFLAGS  -ftree-vectorize "  LDFLAGS="\$LDFLAGS  -ftree-vectorize " ' options];
        if verbose>=0, disp(str); end;
        eval(str);
    end
end
cd(cpth) % return to path

%{ 
	=======================================================================================
	FAIR: Flexible Algorithms for Image Registration, Version 2011
	Copyright (c): Jan Modersitzki
	Maria-Goeppert-Str. 1a, D-23562 Luebeck, Germany
	Email: jan.modersitzki@mic.uni-luebeck.de
	URL:   http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
	=======================================================================================
	No part of this code may be reproduced, stored in a retrieval system,
	translated, transcribed, transmitted, or distributed in any form
	or by any means, means, manual, electric, electronic, electro-magnetic,
	mechanical, chemical, optical, photocopying, recording, or otherwise,
	without the prior explicit written permission of the authors or their
	designated proxies. In no event shall the above copyright notice be
	removed or altered in any way.

	This code is provided "as is", without any warranty of any kind, either
	expressed or implied, including but not limited to, any implied warranty
	of merchantibility or fitness for any purpose. In no event will any party
	who distributed the code be liable for damages or for any claim(s) by
	any other party, including but not limited to, any lost profits, lost
	monies, lost data or data rendered inaccurate, losses sustained by
	third parties, or any other special, incidental or consequential damages
	arrising out of the use or inability to use the program, even if the
	possibility of such damages has been advised against. The entire risk
	as to the quality, the performace, and the fitness of the program for any
	particular purpose lies with the party using the code.
	=======================================================================================
	Any use of this code constitutes acceptance of the terms of the above statements
	=======================================================================================
%}