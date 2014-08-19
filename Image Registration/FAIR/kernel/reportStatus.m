%==============================================================================
% (c) Jan Modersitzki 2011/01/02, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%==============================================================================
% 
% reports current setting
%==============================================================================

caller = dbstack; % identify the name of the calling function
caller = caller(min(length(caller),2)).name;
FAIRmessage(caller)
viewImage('disp');
inter('disp');
distance('disp');
trafo('disp');
regularizer('disp');
fprintf('%% %s\n',char(ones(1,80)*'='));