% (c) Fabian Gigengack 2011/04/13 see FAIR.2 and FAIRcopyright.m.
% http://www.uni-muenster.de/EIMI/
%
% GUI for example files

function EXAMPLES_GUI %#ok<*NASGU>

global glbQuit
glbQuit = false;

% Make GUI
hm = figure(123456789);
set(hm, 'MenuBar', 'none', 'Toolbar','none', ...
        'Name', 'FAIR Examples', ...
        'HandleVisibility', 'callback', ...
        'Color', get(0, 'defaultuicontrolbackgroundcolor'), ...
        'Units', 'pixels', ...
        'Position', [0 0 1130 550], ...
        'CloseRequestFcn', @closefcn);

% Move GUI to center position
movegui(hm, 'center');

pelement = @(parent, style, string, position, value, align, callback) ...
    uicontrol('Parent', parent, 'Style', style, 'String', string, ...
    'Units', 'pixels', 'Position', position, 'Value', value, ...
    'HorizontalAlignment', align, 'Callback', callback);


% Run
pelement(hm, 'pushbutton', 'Run file', [240 520 150 20], 1, 'left', @runFile);

% Open
pelement(hm, 'pushbutton', 'Open file', [390 520 150 20], 1, 'left', @openFile);

% Exit button
pelement(hm, 'pushbutton', 'Exit', [10 520 200 20], 1, 'left', @close);


% Read Example files from examples directory
xmpls = dir([FAIRpath filesep 'examples' filesep '*.m']);

% List
list = {xmpls(:).name};
h_xmpl = pelement(hm, 'listbox', list, [240 10 300 500], 1, 'left', @showHelp);

% Help
h_pnl = uipanel('Parent', hm, 'Title', 'Help', 'Units', 'pixels', ...
    'Position', [560 10 550 530]);
h_help = pelement(h_pnl, 'text', '', [10 10 530 500], 1, 'left', @showHelp);
set(h_help, 'FontName', 'Courier New')


% Dimensions
dims = {'';'1D'; '2D'; '3D'};
cdim = pelement(hm, 'popupmenu', dims, [10 10 200 20], 1, 'left', @update);
pelement(hm,'text', 'Dimensions:', [10 30 200 20], 1, 'left', '');

% Dataset
data = {'';'US';'Hands';'HNSP';'PETCT';'3Dbrain';'MRIhead'};
cdat = pelement(hm, 'popupmenu', data, [10 60 200 20], 1, 'left', @update);
pelement(hm,'text', 'Datasets:', [10 80 200 20], 1, 'left', '');

% Distance Measures
dist = {'';'SSD';'SSDweighted';'MI';'NGF'};
cdis = pelement(hm, 'popupmenu', dist, [10 110 200 20], 1, 'left', @update);
pelement(hm,'text', 'Distance Measures:', [10 130 200 20], 1, 'left', '');

% Regularisation
reg = {'';'Curv';'Elas'};
creg = pelement(hm, 'popupmenu', reg, [10 160 200 20], 1, 'left', @update);
pelement(hm,'text', 'Regularisation:', [10 180 200 20], 1, 'left', '');

% Matrixbased / Matrixfree
mat = {'';'mb';'mf'};
cmat = pelement(hm, 'popupmenu', mat, [10 210 200 20], 1, 'left', @update);
pelement(hm,'text', 'Matrixbased / Matrixfree:', [10 230 200 20], 1, 'left', '');

% Transformation
trf = {'';'affine';'rigid';'rotation';'splineTransformation';'translation'};
ctrf = pelement(hm, 'popupmenu', trf, [10 260 200 20], 1, 'left', @update);
pelement(hm,'text', 'Transformation:', [10 280 200 20], 1, 'left', '');

% Interpolation
int = {'';'linearInter';'linearMatlab';'nnInter';'splineInter'};
cint = pelement(hm, 'popupmenu', int, [10 310 200 20], 1, 'left', @update);
pelement(hm,'text', 'Interpolation:', [10 330 200 20], 1, 'left', '');

% Parametric / Non-Parametric
par = {'';'PIR';'RPIR';'NPIR'};
cpar = pelement(hm, 'popupmenu', par, [10 360 200 20], 1, 'left', @update);
pelement(hm,'text', 'Parametric / Non-Parametric:', [10 380 200 20], 1, 'left', '');

% Multilevel
par = {'';'MLIR';'MLPIR';'MLIRlBFGS'};
cpar = pelement(hm, 'popupmenu', par, [10 410 200 20], 1, 'left', @update);
pelement(hm,'text', 'Multilevel:', [10 430 200 20], 1, 'left', '');

% Search
cedt = pelement(hm, 'edit', '', [10 460 145 20], 1, 'left', @update);
pelement(hm,'pushbutton', 'Search:', [160 460 50 20], 1, 'left', @update);

showHelp;

    function close(~, ~)
        glbQuit = true;
        closefcn
    end

    function closefcn(~, ~)
        if glbQuit, delete(hm); end
    end

    % Update function
    function update(~, ~)
    list = struct('name',{});
    idx = 1;
    for j=1:numel(xmpls)
        glb_add = true;
        
        % Dimensions
        glb_add = glb_add & toadd(xmpls(j).name, 'cdim');
        % Dataset
        glb_add = glb_add & toadd(xmpls(j).name, 'cdat');
        % Distance Measures
        glb_add = glb_add & toadd(xmpls(j).name, 'cdis');
        % Regularisation
        glb_add = glb_add & toadd(xmpls(j).name, 'creg');
        % Matrixbased / Matrixfree
        glb_add = glb_add & toadd(xmpls(j).name, 'cmat');
        % Transformation
        glb_add = glb_add & toadd(xmpls(j).name, 'ctrf');
        % Interpolation
        glb_add = glb_add & toadd(xmpls(j).name, 'cint');
        % Parametric / Non-Parametric
        glb_add = glb_add & toadd(xmpls(j).name, 'cpar');
        % Parametric / Non-Parametric
        glb_add = glb_add & toadd(xmpls(j).name, 'cedt');
        
        if glb_add
            list(idx).name = xmpls(j).name;
            idx = idx + 1;
        end
    end
    set(h_xmpl, 'Value', 1);
    set(h_xmpl, 'String', {list(:).name});
    
    showHelp;
    
        function add = toadd(name, categ)
            add = false;
            str = '';
            eval(['str = get(' categ ', ''String'');']);
            cat = '';
            if strcmp(categ, 'cedt')
                if ~isempty(str), cat = str; end
            else
                eval(['cat = str{get(' categ ', ''Value'')};']);
            end
            if strcmp(cat, '')
                add = true;
            end
            if ~isempty(strfind(lower(name), lower(cat)))
                add = true;
            end
        end
    end

    % Show help
    function showHelp(~, ~)
        filename = getFilename;
        if ~isempty(list), set(h_help, 'String', help(filename)); end
    end

    % Run example file
    function runFile(~, ~)
        filename = getFilename;
        evalin('base', filename(1:end-2));
    end

    % Open example file
    function openFile(~, ~)
        edit(getFilename);
    end

    % Get filename
    function filename = getFilename
        list = get(h_xmpl, 'String');
        if isempty(list), filename = []; return; end
        filename = list{get(h_xmpl, 'Value')};
    end
end

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