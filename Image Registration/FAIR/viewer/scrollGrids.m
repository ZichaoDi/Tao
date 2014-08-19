% (c) Fabian Gigengack 2011/04/13, see FAIRcopyright.m.
% http://www.uni-muenster.de/EIMI/
% http://cvpr.uni-muenster.de/
% 
% function scrollGrids(grid, omega, m, varargin)
% 
% Scrollable viewer for one or multiple 3D grid(s)
%
% Input:
% 	grid        - 3-D grid(s)
% 	omega       - domain(s) of grid(s)
%   m           - size(s) of grid(s)
%   varargin    - variable input arguments
%
% Key-Press Functions:
%   'uparrow'   - show next slice
%   'downarrow' - show previous slice
%   'm'         - show middle slice
%   '1'         - show first dimension
%   '2'         - show second dimension
%   '3'         - show third dimension
%   'i'         - interpolation on/off
%   '*'         - double step size
%   '/'         - half step size
%   ','         - half maximum displayed intensity of vol
%   '.'         - double maximum displayed intensity of vol
%   'c'         - change colormap
%   '+'         - finer spacing of grid
%   '-'         - coarser spacing of grid
%
% Scroll-Wheel Function:
%   show next/previous slice
%
% Example:
%   omega1 = [-10 10 -5 25 -20 20]; m1 = [20 15 20];
%   omega2 = [-15 15 -5  5  -5  5]; m2 = [30  5  5];
%   xc1 = getCellCenteredGrid(omega1, m1);
%   xc2 = getCellCenteredGrid(omega2, m2);
%   xc1 = xc1 + .1 * randn(size(xc1));
%   xc2 = xc2 + .2 * randn(size(xc2));
%   vol = rand(m1);
%   scrollGrids({xc1,xc2}, {omega1,omega2}, {m1,m2}, 'vol', vol, 'color', 'k');

function scrollGrids(grid, omega, m, varargin)

if nargin==0
    runMinimalExample; return;
end

if ~iscell(grid),  grid = {grid}; end
if ~iscell(omega), omega = {omega}; end
if ~iscell(m),     m = {m}; end
nimages = length(grid);
ndomain = length(omega);
nsize   = length(m);
if nimages~=ndomain || nimages~=nsize, error('Different number of volumes, domains and sizes.'); end

dimens  = 3;            % shown dimension
dointer = false;        % interpolate volumes
vol     = [];           % background image (same domain and size as first grid)
mini    = 0;            % minimum value of vol
maxi    = 1;            % minimum value of vol
cmap    = 'jet';        % e.g. 'gray', jet', 'bone', 'hot', 'cool', ...
cbar    = false;        % show colorbar
slice   = round((omega{1}(2*(dimens-1)+2) + omega{1}(2*(dimens-1)+1)) / 2); % shown slice number
stepsz  = 1;            % step size
inter   = @linearInter; % function for image interpolation
spacing = [1 1 1];      % grid spacing
color   = 'b';          % color of grid
minmax  = @max;         % scale to higher (@max) or lower (@min) resolution

for i=1:2:length(varargin) % overwrites default parameter
    eval([varargin{i},'=varargin{',int2str(i+1),'};']);
end

argcheck = @(name) ~isempty(strfind([varargin{1:2:end}],name));
% Check for background image
if argcheck('vol') && ~argcheck('mini'), mini = min(vol(:)); end
if argcheck('vol') && ~argcheck('maxi'), maxi = max(vol(:)); end
if ~argcheck('slice') % Check for predefined slice
    slice = round((omega{1}(2*(dimens-1)+2) + omega{1}(2*(dimens-1)+1)) / 2); % shown slice number
end

% Get combined domain
omegas = omega{1};
for i=2:ndomain, omegas(1:2:end) = min(omegas(1:2:end), omega{i}(1:2:end)); end
for i=2:ndomain, omegas(2:2:end) = max(omegas(2:2:end), omega{i}(2:2:end)); end
% Get combined image size
ms = m{1};
for i=2:nsize, ms = minmax(ms, m{i}); end

if sum(ms~=1)~=3, error('Unknown volume size.'); end
if dimens>numel(ms), dimens = 1; disp('Shown dimension set to 1.'); end
if maxi==mini, maxi = maxi + eps; end

dim = numel(ms);

id = figure;
set(id, 'WindowScrollWheelFcn', @scrollWheelFcn, ...
    'KeyPressFcn', @keyPressFcn, 'CurrentAxes',gca);

doplot

function doplot
    idxO = @(i) (dimens-1)*2+(i);
    idxM = 1:3~=dimens;
    % Get size of axes
    pos = get(get(id, 'CurrentAxes'), 'Position') .* get(id, 'Position');
    pos = pos([3 4]);
    sfac = min(pos ./ ms(idxM));
    clf
    for k=1:nimages
        hold on
        % Extract deformation of k-th grid
        xc = getCellCenteredGrid(omega{k}, m{k});
        d  = reshape(grid{k}(:)-xc, [m{k} dim]);
        slice = min(omegas(idxO(2)), max(omegas(idxO(1)), slice));
        % Get size of grid and displayed slice
        mSlice  = m{k};
        mSliceG = m{k}(idxM);
        mSlice(dimens) = 1;
        % Get domain of grid and displayed slice
        omegaSlice  = omega{k};
        omegaSliceG = omega{k}(1:6~=idxO(1) & 1:6~=idxO(2));
        omegaSlice(idxO(1:2)) = [slice-.5 slice+.5];
        % Get the two dimensions of the deformation that will be displayed
        gridSlice = d(:,:,:,idxM);
        % Get grid for the interpolation of the deformation
        xcSliceG = getCellCenteredGrid(omegaSlice, mSlice);
        % Scale image size
        if dointer || sfac<1, mSlice(idxM) = round(sfac * mSlice(idxM)); end
        % Get grid for the interpolation of the image vol
        xcSlice = getCellCenteredGrid(omegaSlice, mSlice);
        % Display image vol
        if ~isempty(vol) && k==1
            viewImage2D(squeeze(reshape(size(colormap(cmap), 1) * ...
                scale_image(inter(vol, omega{1}, xcSlice), mini, maxi), ...
                mSlice)), omegaSliceG, mSlice(idxM), 'colormap', cmap);
        end
        title(sprintf('Slice %d / %dx%dx%d',slice,ms(1),ms(2),ms(3)));
        if slice>=omega{k}(idxO(1)) && slice<=omega{k}(idxO(2))
            % Inteprolate deformation
            gridSlice = [inter(gridSlice(:,:,:,1), omega{k}, xcSliceG); ...
                         inter(gridSlice(:,:,:,2), omega{k}, xcSliceG)];
            xc = getCellCenteredGrid(omegaSliceG, mSliceG);
            % Plot grid
            plotGrid(xc+gridSlice(:), omegaSliceG,  mSliceG, ...
                'spacing', spacing, 'color', color);
        end
        hold off
    end
    axis off image
    axis ij
    if cbar, colorbar, end
    drawnow
end

function scrollWheelFcn(~,evnt)
    % Scroll wheel events
    slice = slice + stepsz * evnt.VerticalScrollCount;
    doplot;
end

function keyPressFcn(~,evnt)
    % Key press events
    switch evnt.Key
        case 'uparrow'
            slice = slice + stepsz;
        case 'downarrow'
            slice = slice - stepsz;
        case 'multiply'
            stepsz = stepsz * 2;
        case 'divide'
            stepsz = stepsz / 2;
        case '1'
            dimens = 1;
        case '2'
            dimens = 2;
        case '3'
            dimens = 3;
        case 'i'
            dointer = ~dointer;
        case 'm'
            slice = round((omegas(2*(dimens-1)+2) + omegas(2*(dimens-1)+1)) / 2);
        case 'comma'
            maxi = maxi - (maxi - mini) / 2;
        case 'period'
            maxi = maxi + (maxi - mini);
        case 'c'
            switch cmap
                case 'jet'
                    cmap = 'gray';
                case 'gray'
                    cmap = 'hot';
                case 'hot'
                    cmap = 'jet';
            end
        case 'add'
            spacing  = max(1, spacing - 1);
        case 'subtract'
            spacing  = min(ms, spacing + 1);
    end
    doplot
end

function scaled = scale_image(image, mini, maxi)
    % Scale image between 0 and 1 according to mini and maxi
    scaled = (double(image) - mini) / (maxi - mini);
    scaled(scaled<0) = 0;
    scaled(scaled>1) = 1;
end
end

function runMinimalExample
omega1 = [-10 10 -5 25 -20 20]; m1 = [20 15 20];
omega2 = [-15 15 -5  5  -5  5]; m2 = [30  5  5];
xc1 = getCellCenteredGrid(omega1, m1);
xc2 = getCellCenteredGrid(omega2, m2);
xc1 = xc1 + .1 * randn(size(xc1));
xc2 = xc2 + .2 * randn(size(xc2));
vol = rand(m1);
scrollGrids({xc1,xc2}, {omega1,omega2}, {m1,m2}, 'vol', vol, 'color', 'k');
help(mfilename)
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