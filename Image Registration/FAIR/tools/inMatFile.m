% -----------------------------------------------------------------------------------
% (c) Jan Modersitzki 2011-01-11, see FAIR.2 and FAIRcopyright.m.
% checks whether var is contained in matfile
% -----------------------------------------------------------------------------------

function OK = inMat(matfile,var)
what = whos('-file',matfile);
OK = 1;
for k=1:length(var),
  if isempty(find(strcmp({what(:).name},var{k})==1));
    OK = 0; error(var{k}); 
  end;
end;