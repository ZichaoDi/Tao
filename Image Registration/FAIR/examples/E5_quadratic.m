%==============================================================================
% (c) Jan Modersitzki 2010/12/28, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% Tutorial for FAIR:  Quadratic LM registration for hand example
%
% - load data (see setupHandData)
% - compute least squares solution and plot
% - this file for demonstration, see also LMreg for a convenient variant
%==============================================================================

clear, close all, help(mfilename)

setupHandData; xc = reshape(getCellCenteredGrid(omega,m),[],2); %LM(7,:) = [];
Q = [ones(size(LM,1),1),LM(:,[3:4]),LM(:,3).^2,LM(:,4).^2,LM(:,3).*LM(:,4)];
wc = (Q'*Q)\(Q'*LM(:,1:2));

quad = @(w,x1,x2) (w(1)+w(2)*x1+w(3)*x2+w(4)*x1.^2+w(5)*x2.^2+w(6)*x1.*x2); 
yc = [quad(wc(:,1),xc(:,1),xc(:,2));quad(wc(:,2),xc(:,1),xc(:,2))];
LM(:,[5,6]) = [quad(wc(:,1),LM(:,3),LM(:,4)),quad(wc(:,2),LM(:,3),LM(:,4))];

P5_LM; % for nice plots
