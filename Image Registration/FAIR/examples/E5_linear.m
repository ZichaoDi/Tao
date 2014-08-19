%==============================================================================
% (c) Jan Modersitzki 2010/12/28, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% Tutorial for FAIR:  Linear LM registration for hand example
%
% - load data (see setupHandData)
% - compute least squares solution and plot
% - this file for demonstration, see also LMreg for a convenient variant
%==============================================================================

clear, close all, help(mfilename)

setupHandData; xc = reshape(getCellCenteredGrid(omega,m),[],2);
Q  = [LM(:,3:4),ones(size(LM,1),1)];
wc = (Q'*Q)\(Q'*LM(:,1:2));
yc  = [(wc(1,1)*xc(:,1) + wc(2,1)*xc(:,2) + wc(3,1)),...
       (wc(1,2)*xc(:,1) + wc(2,2)*xc(:,2) + wc(3,2)) ];
LM(:,[5,6]) = ...    
     [(wc(1,1)*LM(:,3) + wc(2,1)*LM(:,4) + wc(3,1)),...
      (wc(1,2)*LM(:,3) + wc(2,2)*LM(:,4) + wc(3,2)) ];

P5_LM; % for nice plots

