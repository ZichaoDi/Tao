% (c) Lars Ruthotto 2011/02/08, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de
% 
% function scrollVolume(T,omega,m,varargin)
%
% scroll and compare 3D images
%
% Handling:
% ---------
% Mousewheel (up/down) - scroll to next slices
% Arrow      (up/down) - scroll to next slices
% Slider     (up/down) - scroll to next slices
% Choose view and colormap by GUI
%
% Example:
% --------
% >> load mice3D; scrollVolume(dataR,omega,m,dataT);
%
% Input:
% ------
%   T       - 3D image 
%   omega   - representation of computational domain
%   m       - discretization size
%   varagin - additional images for comparison

function scrollVolume(T,omega,m,varargin)

if nargin==0
    runMinimalExample; return;
end

view = 'xy'; % scroll z-axis
if ~exist('m',    'var'), m     = size(T); end;
if ~exist('omega','var'), 
  omega = reshape([ones(1,length(m));m],1,[]); 
end;
% load images
nImg = length(varargin)+size(T,4);
for k=1:size(T,4),
    im{k} = T(:,:,:,k);
end
for k=size(T,4)+1:nImg,% include additional images
    im{k} = varargin{k-1};
    if any(size(im{k}) ~= m),
        error('All image must have same dimensions');
    end;
end
% get image names
imname  = cell(nImg,1);
for k=1:size(T,4),
    imname{k} = inputname(1);
end

for k=size(T,4)+1:nImg
    imname{k} = inputname(k+2);
end
% set initial slices
slice = round(m(3)/2);

p = round(sqrt(nImg)); % number of rows
q = ceil(nImg/p);      % number of columns
if (p-1)*q >= nImg, p=p-1; end


fig = figure;
clf;
set(fig,  'WindowScrollWheelFcn',@figMouseWheelScroll, ...
          'KeyPressFcn',@figUpDownKeyScroll,...
          'WindowButtonMotionFcn',@sliderCallback);

him  = zeros(nImg,1);
for k =1:nImg,
    R = im{k};
    subplot(p,q,k);
    [I,om,mc] = getSlice(R,slice);
    him(k)  = viewImage2Dsc(I,om,mc);
    colormap gray;
    imgTitle(k);
end

hslider=uicontrol( ...
    'Style','slider', ...
    'Units','normalized', ...
    'position',[0.99 0 0.2 1 ], ...
    'value',slice, ...
    'min',1, ...
    'max',getNumSlices);
set(hslider,'callback',@sliderCallback);

% UIcontrol to select view={'xy','xz','yz'}
hview = uicontrol('Style', 'popup',...
       'String', 'xy|xz|yz',...
       'Position', [0 0 100 30]);
set(hview,'callback',@viewCallback);

%UIcontrol to select colormap={'gray','jet','gray(inv)'}
hcmap = uicontrol('Style', 'popup',...
       'String', 'gray|jet|gray(inv)',...
       'Position', [100 0 100 30]);
set(hcmap,'callback',@cmapCallback);
drawnow

    function cmapCallback(hcmap,eventdata)
        str = get(hcmap,'string');
        map_idx = get(hcmap,'Value');
        map = str(map_idx,:);
        bInv = 0;
        if findstr(map,'(inv)')
            map = regexprep(map,'(inv)','');
            bInv = 1;
        end
        C = eval(['colormap(' map ')']);
        if bInv,
            C = flipud(C);
            colormap(C);
        end
        drawnow;
    end
    function viewCallback(hpop,eventdata)
        switch get(hpop,'Value')
            case 1
                view = 'xy';
            case 2
                view = 'xz';
            case 3
                view = 'yz';
        end
            drawSlice;
    end
    function sliderCallback(src,event)
        slice = round(get(hslider,'value'));
        drawSlice;
    end

    function figMouseWheelScroll(src,evnt)
        if evnt.VerticalScrollCount > 0
            if slice<getNumSlices, slice = slice+1; drawSliceAndUpdateSlider; end
        elseif evnt.VerticalScrollCount < 0
            if slice>1,  slice = slice-1;  drawSliceAndUpdateSlider;  end
        end
    end

    function figUpDownKeyScroll(src,evnt)
        if strcmp(evnt.Key,'uparrow')
            if slice<getNumSlices, slice = slice+1; drawSliceAndUpdateSlider; end
        elseif strcmp(evnt.Key,'downarrow')
            if slice>1,  slice = slice-1;   drawSliceAndUpdateSlider;  end
        end
    end

    function drawSliceAndUpdateSlider
        set(hslider,'value',slice);
        drawSlice;
    end
    function [slice,om,mc,ns] = getSlice(R,slice)
        switch view
            case 'xy'
                slice = squeeze(R(:,:,slice));
                mc = m(1:2);
                om = omega(1:4);
            case 'xz'
                slice = squeeze(R(:,slice,:));
                mc = m([1 3]);
                om = omega([1 2 5 6]);
            case 'yz'
                slice = squeeze(R(slice,:,:));
                mc = m(2:3);
                om = omega(3:6);
        end 
    end
    function ns = getNumSlices
        switch view
            case 'xy'
                ns = m(3);
            case 'xz'
                ns = m(2);
            case 'yz'
                ns = m(1);
        end 
    end
    function imgTitle(k)
        switch view
            case 'xy'
                str = 'z';
            case 'xz'
                str = 'y';
            case 'yz'
                str = 'x';
        end
        title(sprintf('%s, %s-slice %d of %d [%d x %d]',imname{k},str,slice,...
                       getNumSlices,mc(1),mc(2)))
    end
    function drawSlice
        slice = min(getNumSlices,slice);
        for k =1:nImg,
            R = im{k};
            [I, om,mc] = getSlice(R,slice);
            ax(k) = subplot(p,q,k);
            him(k)  = viewImage2Dsc(I,om,mc);
            imgTitle(k)
        end
        drawnow;
    end
end

function runMinimalExample
load mice3D; scrollVolume(dataR,omega,m,dataT);
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