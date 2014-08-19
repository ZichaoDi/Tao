%==============================================================================
% (c) Jan Modersitzki 2011/01/02, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% Tutorial for FAIR:  SSD versus rotation
%
%   - data                 PETCT, Omega=(0,140)x(0,151), level=6, m=[128,128]
%   - viewer               viewImage2D
%   - interpolation        splineInter
%   - distance             SSD
%   - transformation       rotation2D
%   - optimization         Gauss-Newton
%==============================================================================

clear, close all, help(mfilename);

setupPETCTdata; level = 6; omega = MLdata{level}.omega; m = MLdata{level}.m;
viewImage('reset',viewOptn{:},'axis','off');
inter('reset','inter','splineInter','regularizer','moments','theta',1e0);
[T,R] = inter('coefficients',MLdata{level}.T,MLdata{level}.R,omega,'out',0);
distance('reset','distance','SSD');
center = (omega(2:2:end)-omega(1:2:end))'/2;
trafo('reset','trafo','rotation2D','c',center);
fprintf('%20s : %s\n','viewImage',viewImage);
fprintf('%20s : %s\n','inter',inter);
fprintf('%20s : %s\n','distance',distance);
fprintf('%20s : %s\n','trafo',trafo);

xc = getCellCenteredGrid(omega,m);
Rc = inter(R,omega,xc);
wc = linspace(-pi/2,pi/2,51);
Dc = zeros(size(wc));

% run the loop over all rotations
for j = 1:length(wc),
  yc = trafo(wc(j),xc);                 % compute transformed grid
  Tc = inter(T,omega,yc);               % compute transformed image
  [Dc(j),rc] = distance(Tc,Rc,omega,m); % compute distance
  
  % visualize
  if j == 1,
    FAIRfigure(1,'figname',mfilename); clf;
    subplot(2,2,1);  viewImage(Rc,omega,m);            th(1) = title('R');
    subplot(2,2,2);  vh = viewImage(Tc,omega,m);       th(2) = title('T(Y)');
    subplot(2,2,3);  rh = viewImage(128+rc/2,omega,m); th(3) = title('|T(Y)-R|');
    subplot(2,2,4);  ph = plot(wc(1),Dc(1),'r.','markersize',20);
    th(4) = title(sprintf('%s versus rotation',distance));
    axis([wc(1),wc(end),-inf,inf]); hold on;
    set(th,'fontsize',30);
    pause;
  else
    set(vh,'cdata',reshape(Tc,m)')
    set(rh,'cdata',reshape(128+rc/2,m)');
    subplot(2,2,4); set(ph,'visible','off');
    plot(wc(1:j),Dc(1:j),'k-','linewidth',2);
    ph = plot(wc(j),Dc(j),'r.','markersize',20);
    drawnow, pause(1/100)
  end;
  fprintf('.'); if ~rem(j,50) || j == length(wc), fprintf('\n'); end;
end;
