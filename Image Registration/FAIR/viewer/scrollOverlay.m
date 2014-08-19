% (c) Fabian Gigengack 2011/04/13, see FAIRcopyright.m.
% http://www.uni-muenster.de/EIMI/
% http://cvpr.uni-muenster.de/
% 
% function scrollOverlay(vol1, omega1, vol2, omega2, varargin)
% 
% Scrollable overlay function for two 3D volumes
%
% Input:
% 	vol1        - 3-D volume
% 	omega1      - domain of vol1
% 	vol2        - 3-D volume
% 	omega2      - domain of vol2
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
%   'a'         - switch to amide view (amide.sourceforge.net)
%   '*'         - double step size
%   '/'         - half step size
%   ','         - half maximum displayed intensity of vol1
%   '.'         - double maximum displayed intensity of vol1
%   'k'         - half maximum displayed intensity of vol2
%   'l'         - double maximum displayed intensity of vol2
%   '-'         - fade to vol1
%   '+'         - fade to vol2
%
% Scroll-Wheel Function:
%   show next/previous slice
%
% Example1:
% 	load mice3D;
%   scrollOverlay(dataR, omega, dataT, omega, 'dointer', true, ...
%                 'mini1', 0, 'maxi1', 200, 'mini2', 0, 'maxi2', 200);
%
% Example2:
%   load brain3D;
% 	scrollOverlay(dataR, omega, dataT, omega, 'stepsz',  20/128, ...
%                 'dointer', true, 'cmap1', 'gray', 'cmap2', 'gray');

function scrollOverlay(vol1, omega1, vol2, omega2, varargin)

if nargin==0
    runMinimalExample; return;
end

% The following arguments can be overloaded with varargin
dimens  = 3;            % shown dimension
doamide = false;        % amide view (amide.sourceforge.net)
dointer = false;        % interpolate volumes
mini1   = min(vol1(:)); % minimum value of vol1
maxi1   = max(vol1(:)); % maximum value of vol1
mini2   = min(vol2(:)); % minimum value of vol2
maxi2   = max(vol2(:)); % maximum value of vol2
show    = @imagesc;     % @imshow or @imagesc
cmap1   = 'jet';        % e.g. 'gray', jet', 'bone', 'hot', 'cool', ...
cmap2   = 'jet';        % e.g. 'gray', jet', 'bone', 'hot', 'cool', ...
cbar    = false;        % show colorbar
fac     = .5;           % blending factor
slice   = round((omega1(2*(dimens-1)+2) + omega1(2*(dimens-1)+1)) / 2); % shown slice number
stepsz  = 1;            % step size
inter   = @linearInter; % function for image interpolation

for i=1:2:length(varargin) % overwrites default parameter
	eval([varargin{i},'=varargin{',int2str(i+1),'};']);
end

argcheck = @(name) isempty(strfind([varargin{1:2:end}],name));
if argcheck('slice') % Check for predefined slice
    slice = round((omega1(2*(dimens-1)+2) + omega1(2*(dimens-1)+1)) / 2); % shown slice number
end

m1 = size(vol1);
m2 = size(vol2);

if sum(m1~=1)~=3 || sum(m2~=1)~=3, error('Unknown volume size.'); end
if dimens>numel(m1), dimens = 1; disp('Shown dimension set to 1.'); end
if maxi1==mini1, maxi1 = maxi1 + eps; end
if maxi2==mini2, maxi2 = maxi2 + eps; end

id = figure;
set(id, 'WindowScrollWheelFcn', @scrollWheelFcn, ...
    'KeyPressFcn', @keyPressFcn, 'CurrentAxes',gca);

hh = (omega1(2:2:end)-omega1(1:2:end))./size(vol1);
hh = min(hh, (omega2(2:2:end)-omega2(1:2:end))./size(vol2));

doplot

function doplot
    % Get combined domain of vol1 and vol2
    omega = zeros(size(omega1));
    omega(1:2:end) = min(omega1(1:2:end), omega2(1:2:end));
    omega(2:2:end) = max(omega1(2:2:end), omega2(2:2:end));
    % Get combined image size of vol1 and vol2
    m = ceil((omega(2:2:end)-omega(1:2:end))./hh);
    % Get size of axes
    pos = get(get(id, 'CurrentAxes'), 'Position') .* get(id, 'Position');
    if doamide, pos = pos([3 4]); else pos = pos([4 3]); end
    slice = min(omega((dimens-1)*2+2), max(omega((dimens-1)*2+1), slice));
    switch dimens
        case 1
            omega = [slice-.5 slice+.5 omega(3:6)];
            sfac  = min(pos ./ m([2 3]));
        case 2
            omega = [omega(1:2) slice-.5 slice+.5 omega(5:6)];
            sfac  = min(pos ./ m([1 3]));
        case 3
            omega = [omega(1:4) slice-.5 slice+.5];
            sfac  = min(pos ./ m([1 2]));
    end
    % Scale image size in case of interpolation
    if dointer || sfac<1, m = round(sfac * m); end
    m(dimens) = 1;
    % Interpolate vol1 and vol2 (one slice)
    grid = getCellCenteredGrid(omega, m);
    img1 = amide(squeeze(reshape(inter(squeeze(vol1), omega1, grid), m)));
    img2 = amide(squeeze(reshape(inter(squeeze(vol2), omega2, grid), m)));
    % Compute domain of slice
    omega2D = omega([1 1 2 2 3 3]~=dimens);
    if doamide, omega2D = omega2D([3 4 1 2]); end
    m = m(m~=1); if doamide, m = fliplr(m); end
    h2D = (omega2D(2:2:end)-omega2D(1:2:end))./m;
    xi = @(i) (omega2D(2*i-1)+h2D(i)/2:h2D(i):omega2D(2*i)-h2D(i)/2)';
    % Compute overlay and display
    show(xi(2),xi(1),imoverlay(scale_image(img1, mini1, maxi1), ...
                               scale_image(img2, mini2, maxi2), ...
                               fac, eval(cmap1), eval(cmap2)), [0 255]);
    axis image; axis off;
    title(['Slice ' num2str(slice) ' (' num2str(1-fac) ' * im1 + ' num2str(fac) ' * im2)'])
    if cbar, colorbar, end
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
        case 'a'
            doamide = ~doamide;
        case 'i'
            dointer = ~dointer;
        case 'add'
            fac = min(1, fac + .1);
        case 'subtract'
            fac = max(0, fac - .1);
        case 'm'
            slice = round((omega1(2*(dimens-1)+2) + omega1(2*(dimens-1)+1)) / 2);
        case 'comma'
            maxi1 = maxi1 - (maxi1 - mini1) / 2;
        case 'period'
            maxi1 = maxi1 + (maxi1 - mini1);
        case 'k'
            maxi2 = maxi2 - (maxi2 - mini2) / 2;
        case 'l'
            maxi2 = maxi2 + (maxi2 - mini2);
    end
    doplot;
end

function scaled = scale_image(image, mini, maxi)
    % Scale image between 0 and 1 according to mini and maxi
    scaled = (double(image) - mini) / (maxi - mini);
    scaled(scaled<0) = 0;
    scaled(scaled>1) = 1;
end

function out = imoverlay(im1, im2, factor, map1, map2)
    % Overlay im1 with colormap map1 and im2 with colormap map2
    if numel(factor)~=1 || factor<0 || factor>1, factor = 0.5; end
    im1 = double(255 * ind2rgb(uint8(im1*size(map1,1)), map1));
    im2 = double(255 * ind2rgb(uint8(im2*size(map2,1)), map2));
    out = uint8((1 - factor) * im1 + factor * im2);
end

function img = amide(img)
    if ~doamide, return, end
    % switch to amide view (amide.sourceforge.net)
    switch dimens
        case 1
            img = fliplr(rot90(img, -1));
        case 2
            img = fliplr(rot90(img, -1));
        case 3
            img = rot90(img);
    end
end
end

function runMinimalExample
load brain3D;
scrollOverlay(dataR, omega, dataT, omega, 'stepsz',  20/128, ...
    'cmap1', 'gray', 'cmap2', 'gray');
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