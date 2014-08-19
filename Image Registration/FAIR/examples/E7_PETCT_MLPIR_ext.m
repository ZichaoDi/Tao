%==============================================================================
% (c) Jan Modersitzki 2011/01/02, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% Tutorial for FAIR: various distances and Multi-Level Parametric Image Registration
%
%   - data                 PETCT, Omega=(0,140)x(0,151), level=4:7, m=[128,128]
%   - viewer               viewImage2D
%   - interpolation        splineInter
%   - distance             {'SSD','NCC','MI','NGF'}
%   - transformation       affine2D
% see also E7_PETCT_MLPIR
%==============================================================================

clear, close all, help(mfilename);

setupPETCTdata;

FAIRfigure(12,'position','default');
FAIRfigure(13,'position','default');
FAIRfigure(14,'position','default');

FAIRprint('reset','folder',fullfile(FAIRpath,'../temp','PETCT'),...
  'pause',1,'obj','gca','format','jpg','draft','off');
FAIRprint('disp');

mname = mfilename;
File  = @(str) sprintf('%s-%s',mname,str);

outpath = fullfile(FAIRpath,'temp','PETCT');
if ~exist(outpath,'dir'),
  mkdir(outpath);
end;

Write = @(T,str) imwrite(uint8(flipud(reshape(T,m)')),...
  fullfile(outpath,[File(str),'.jpg']));

DM = {'SSD','NCC','MIcc','NGFdot'};
map = @(T) 160*(T/160).^(0.5);

for dm = 1:length(DM),

  inter('reset','inter','splineInter','regularizer','moments','theta',1e0);
  level = 5; omega = MLdata{level}.omega; m = MLdata{level}.m;
  [T,R] = inter('coefficients',MLdata{level}.T,MLdata{level}.R,omega);
  trafo('reset','trafo','affine2D');  % initialize transformation
  wStop = trafo('w0');
  wOpt = wStop;

  distance('reset','distance',DM{dm});
  distance('disp');
  
  if 1,  
    wSmooth =  MLPIR(MLdata,'minLevel',5,'plotIter',0,'plotMLiter',0);

    inter('reset','inter','splineInter','regularizer','moments','theta',1e-3);
    level = length(MLdata); omega = MLdata{level}.omega; m = MLdata{level}.m;
    [T,R] = inter('coefficients',MLdata{level}.T,MLdata{level}.R,omega);
    xc = getCellCenteredGrid(omega,m);
    Rc = inter(R,omega,xc);
    fctn = @(wc) PIRobjFctn(T,Rc,omega,m,0,[],[],xc,wc);
    wOpt = GaussNewton(fctn,wSmooth,'yStop',wStop);
  end;

  omega = MLdata{end}.omega; m = MLdata{end}.m;
  xc = getCellCenteredGrid(omega,m);
  R0 = inter(R,omega,xc);
  T0 = inter(T,omega,xc);

  if dm == 0,
    FAIRfigure(13)
    str = 'R0'; var = R0;
    viewImage(var,omega,m,'axis','off','title',File(str));  Write(var,str);

    FAIRfigure(12)
    str = 'T0'; var = T0;
    viewImage(var,omega,m,'axis','off','title',File(str));  Write(var,str);

    FAIRfigure(14); clf;
    Tk = rescale(T0,0,255,3);
    overlayImage2D(Tk,R0,omega,m);
    FAIRprint(File(['D0']));
    return
  end;

  Yopt = trafo(wOpt,xc);
  Tc = inter(T,omega,Yopt);

  FAIRfigure(12); clf;
  str = ['G-',distance];  var = T0;
  viewImage(var,omega,m,'axis','off'); hold on;
  plotGrid(grid2grid(Yopt,m,'centered','nodal'),omega,m,...
    'spacing',ceil(m/16),'linewidth',3,'color','w');
  set(gca,'position',[0 0 1 1]);
  FAIRprint(File(['G-',distance]));

  FAIRfigure(13); clf;
  str = ['T-',distance]; var = Tc;
  viewImage(var,omega,m,'axis','off','title',File(str));  Write(var,str);

  FAIRfigure(14); clf;
  Tk = rescale(Tc,0,255,3);
  overlayImage2D(Tk,R0,omega,m);
  FAIRprint(File(['D-',distance]));
end;

pause(5); close all
