% (c) Fabian Gigengack 2011/04/13, see FAIRcopyright.m.
% http://www.uni-muenster.de/EIMI/
% http://cvpr.uni-muenster.de/
% 
% function scrollImages(vol, omega, varargin)
% 
% Scrollable viewer for one or multiple 3D volume(s)
%
% Input:
% 	vol         - 3-D volume(s)
% 	omega       - domain(s) of vol
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
%   ','         - half maximum displayed intensity of vol
%   '.'         - double maximum displayed intensity of vol
%
% Scroll-Wheel Function:
%   show next/previous slice
%
% Example1:
%   load brain3D;
%   scrollImages({dataT, dataR}, {omega, omega}, 'stepsz',  20/128, 'cmap', 'gray');
%
% Example2:
%   load brain3D;
%   scrollImages(dataT, omega, 'stepsz',  20/128, 'cmap', 'gray');

function scrollImages(vol, omega, varargin)

if nargin==0
    runMinimalExample; return;
end

if ~iscell(vol), vol = {vol}; end
if ~iscell(omega), omega = {omega}; end
nimages = length(vol);
ndomain = length(omega);
if nimages~=ndomain, error('Different number of volumes and domains.'); end

% The following arguments can be overloaded with varargin
dimens  = 3;            % shown dimension
doamide = false;        % amide view (amide.sourceforge.net)
dointer = false;        % interpolate volumes
mini    = min2(vol{1}); % minimum value of vol
for i=2:nimages, mini = min(mini, min2(vol{i})); end
maxi    = max2(vol{1}); % maximum value of vol
for i=2:nimages, maxi = max(maxi, max2(vol{i})); end
show    = @imagesc;     % @imshow or @imagesc
cmap    = 'jet';        % e.g. 'gray', jet', 'bone', 'hot', 'cool', ...
cbar    = false;        % show colorbar
slice   = round((omega{1}(2*(dimens-1)+2) + omega{1}(2*(dimens-1)+1)) / 2); % shown slice number
stepsz  = 1;            % step size
inter   = @linearInter; % function for image interpolation
row     = -1;           % number of rows
col     = -1;           % number of columns
minmax  = @min;         % scale to higher (@max) or lower (@min) resolution

for i=1:2:length(varargin) % overwrites default parameter
	eval([varargin{i},'=varargin{',int2str(i+1),'};']);
end

argcheck = @(name) isempty(strfind([varargin{1:2:end}],name));
if argcheck('slice') % Check for predefined slice
    slice = round((omega{1}(2*(dimens-1)+2) + omega{1}(2*(dimens-1)+1)) / 2); % shown slice number
end

m = size(vol{1});

if sum(m~=1)~=3, error('Unknown volume size.'); end
if dimens>numel(m), dimens = 1; disp('Shown dimension set to 1.'); end
if maxi==mini, maxi = maxi + eps; end

p = round(sqrt(nimages));
q = ceil(nimages/p);
if col~=-1 && mod(nimages,col)==0, q=col; p=nimages/col; end
if row~=-1 && mod(nimages,row)==0, q=nimages/row; p=row; end
if (p-1)*q >= nimages, p=p-1; end

id = figure;
set(id, 'WindowScrollWheelFcn', @scrollWheelFcn, ...
    'KeyPressFcn', @keyPressFcn, 'CurrentAxes',gca);

% Get combined domain
omegas = omega{1};
for i=2:ndomain, omegas(1:2:end) = min(omegas(1:2:end), omega{i}(1:2:end)); end
for i=2:ndomain, omegas(2:2:end) = max(omegas(2:2:end), omega{i}(2:2:end)); end
% Get combined spacing
h  = (omega{1}(2:2:end)-omega{1}(1:2:end))./size(vol{1});
for i=2:nimages, h  = minmax(h, (omega{i}(2:2:end)-omega{i}(1:2:end))./size(vol{i})); end

doplot

function doplot
    omegaSlice = omegas;
    % Get size of axes
    pos = get(get(id, 'CurrentAxes'), 'Position') .* get(id, 'Position');
    if doamide, pos = pos([3 4]); else pos = pos([4 3]); end
    slice = min(omegaSlice((dimens-1)*2+2), max(omegaSlice((dimens-1)*2+1), slice));
    mSlice = ceil((omegas(2:2:end)-omegas(1:2:end))./h);
    switch dimens
        case 1
            omegaSlice = [slice-.5 slice+.5 omegaSlice(3:6)];
            sfac  = min(pos ./ mSlice([2 3]));
        case 2
            omegaSlice = [omegaSlice(1:2) slice-.5 slice+.5 omegaSlice(5:6)];
            sfac  = min(pos ./ mSlice([1 3]));
        case 3
            omegaSlice = [omegaSlice(1:4) slice-.5 slice+.5];
            sfac  = min(pos ./ mSlice([1 2]));
    end
    % Scale image size in case of interpolation
    if dointer || sfac<1, mSlice = round(sfac * mSlice); end
    mSlice(dimens) = 1;
    % Interpolate vol1 and vol2 (one slice)
    grid = getCellCenteredGrid(omegaSlice, mSlice);
    mm = mSlice(mSlice~=1); if doamide, mm = fliplr(mm); end
    % Compute overlay
    vol_tmp = zeros(mm(1)*p, mm(2)*q);
    for k=1:nimages
        [cc, rr] = ind2sub([q p], k);
        vol_tmp(mm(1)*(rr-1)+(1:mm(1)),mm(2)*(cc-1)+(1:mm(2))) = ...
            amide(squeeze(reshape(inter(vol{k}, omega{k}, grid), mSlice)));
    end
    % Compute domain of slice
    omega2D = omegas([1 1 2 2 3 3]~=dimens);
    if doamide, omega2D = omega2D([3 4 1 2]); end
    omegaTotal = omega2D .* [p p q q];
    h2D = (omega2D(2:2:end)-omega2D(1:2:end))./mm;
    xi = @(i) (omegaTotal(2*i-1)+h2D(i)/2:h2D(i):omegaTotal(2*i)-h2D(i)/2)';
    % Display overlay
    show(xi(2),xi(1),vol_tmp, [mini maxi]);
    colormap(cmap); axis image; axis off;
    title(['Slice ' num2str(slice)])
    if cbar, colorbar, end
    
    % Plot grid between images
    hold on
    clr = get(id, 'Color');
    for pp=1:p-1
        plot([omegaTotal(3) omegaTotal(4)], pp*[omega2D(2) omega2D(2)], 'Color', clr)
    end
    for qq=1:q-1
        plot(qq*[omega2D(4) omega2D(4)], [omegaTotal(1) omegaTotal(2)], 'Color', clr)
    end
    hold off
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
        case 'm'
            slice = round((omegas(2*(dimens-1)+2) + omegas(2*(dimens-1)+1)) / 2);
        case 'comma'
            maxi = maxi - (maxi - mini) / 2;
        case 'period'
            maxi = maxi + (maxi - mini);
    end
    doplot;
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

function out = min2(in)
    out = min(in(:));
end

function out = max2(in)
    out = max(in(:));
end
end

function runMinimalExample
load brain3D;
scrollImages({dataT, dataR}, {omega, omega}, 'stepsz',  20/128, 'cmap', 'gray');
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