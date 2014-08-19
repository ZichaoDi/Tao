%==============================================================================
% (c) Jan Modersitzki 2010/12/27, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% Tutorial for FAIR: 3D data ordering
%
%==============================================================================

Tij=reshape(1:12,[3,2,2]), 
disp(flipdim(flipdim(permute(Tij,[3,2,1]),3),1))
disp(['T(:)''= [ 9  3 12  6  8  2 11  5  7  1 10  4]'])
