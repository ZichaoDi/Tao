%==============================================================================
% (c) Jan Modersitzki 2010/12/27, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% Tutorial for FAIR: Used to create Fig. 3.13 p.40
%
%==============================================================================

setupUSData;
omega = MLdata{end}.omega; m = MLdata{end}.m;
distance('reset','distance','SSD');
inter('reset','inter','splineInter');
lMax = length(MLdata);
m  = @(l) MLdata{l}.m;
xc = @(l) getCellCenteredGrid(omega,m(l));

Name  = @(c,l) sprintf('ML%c-%d',c,l);
Write = @(c,l,R,m) imwrite(uint8(flipud(reshape(R,m)')),...
  fullfile(FAIRpath,'temp',[Name(c,l),'.jpg']));

for level=3:length(MLdata);
  R = inter('coefficients',MLdata{level}.T,[],omega);
  Rc = inter(R,omega,xc(level));
  Rf = inter(R,omega,xc(lMax));
  figure(level); clf;
  subplot(1,2,1); viewImage(Rf,omega,m(lMax));
  subplot(1,2,2); viewImage(Rc,omega,m(level)); pause(1/6);
  Write('f',level,Rf,m(lMax));
  Write('c',level,Rc,m(level));
end;
