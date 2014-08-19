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
% see also E7_PETCT_MLPIR_ext
%==============================================================================

% =========================================================================
function E7(example);

if nargin == 0, example = 'PETCT';  end;
imagePath  = fullfile(FAIRpath,'temp',example);
switch example
  case 'HNSP'
    setupHNSPdata;
    theta = 1;
  case 'MRIhead',
    setupMRIdata;
    theta = 1;
  case 'PETCT',
    setupPETCTdata;
    theta = 0;
  otherwise, return;
end
level = 6; omega = MLdata{level}.omega; m = MLdata{level}.m;

DM = {'SSD','NCC','MIcc','NGFdot'};
viewImage('set',viewOptn{:})
filename = fullfile(FAIRpath,'temp',sprintf('%s-%s-%s.mat',mfilename,'rotation',example));

% 2e1 MRIhead
inter('reset','inter','splineInter','regularizer','moments','theta',theta);
[T,R] = inter('coefficients',MLdata{level}.T,MLdata{level}.R,omega);
xc    = getCellCenteredGrid(omega,m); 
xc    = reshape(xc,[],2);
Tin   = inter(T,omega,xc);
Rin   = inter(R,omega,xc);

%%  plots
k = figure(1); clf; set(k,'position',FAIRposition(800),'color','w');
subplot(1,2,1); viewImage(Tin,omega,m);
subplot(1,2,2); viewImage(Rin,omega,m);
pause(1);

center = (omega(2:2:end)-omega(1:2:end))'/2;
trafo('reset','trafo','rotation2D','c',center);
trafo('w0');
if ~exist(filename,'file'), 
  w = pi/4*linspace(-1,1,101)';
  save(filename,'w','DM');  
else
  clear DM*
  load(filename)
end;

edge = 100;

for k=1:length(DM),
  variable = ['DM',DM{k}];
  var = whos('-file',filename);
  j = find(strcmp({var(:).name},variable)==1);

  if isempty(j),
    fprintf('============== %s ====================\n\n',variable)
    disp([variable,'=dm;']);
    dm = zeros(size(w));
    for j=1:length(w),
      Y  = trafo(w(j),xc(:));
      Tc = inter(T,omega,Y);
      dm(j) = feval(DM{k},Tc,Rin,omega,m,'edge',edge);
      if j== 1,
        figure(3);
        ph = viewImage(Tc,omega,m);
      else
        set(ph,'cdata',reshape(Tc,m)'); drawnow
        pause(1/100)
      end;
      title(sprintf('%d/%d:%d/%d',k,length(DM),j,numel(w)));
    end;
    eval([variable,'=dm;']);
    save(filename,'-append',variable);
  end;
end;

%%
load(filename)


filename = fullfile(FAIRpath,'temp',...
  sprintf('%s-%s-%s.mat',mfilename,'rotation',example));

FAIRprint('set','folder',imagePath,'prefix',[example,'-rotation-'],...
  'pause',1,'obj','gcf','format','jpg');
FAIRprint('disp');

rotPrint = @(k) FAIRprint(sprintf('%s',DM{k}));

for k=1:4,
  variable = ['DM',DM{k}];
  eval(['dm=',variable,';']);
  dm = dm(1:length(w));
  [wOpt,j] = min(dm);
  figure(k); clf; set(k,'position',FAIRposition(800),'color','w');
  ph = plot(w,dm,'-',w(j),dm(j),'*');
  set(ph,'linewidth',2,'color','k','markersize',20);
  set(gca,'fontsize',30);
  a = max(dm)-0.2*(max(dm)-min(dm));
  th = text(w(j),a,['$w^*=',num2str(w(j)),'$']);
  set(th,'fontsize',30,'interpreter','latex','horizontalalignment','center');
  rotPrint(k);
end;


%% finish
filename = fullfile(FAIRpath,'temp',...
  sprintf('%s-%s-%s.mat',mfilename,'translation',example));

FAIRprint('set','folder',imagePath,'prefix',[example,'-translation-'],...
  'pause',1,'obj','gcf','format','jpg');
FAIRprint('disp');

trafo('reset','trafo','translation2D');
if ~exist(filename,'file'), 
  [w1,w2] = ndgrid(0.1*(omega(2)-omega(1))*linspace(-1,1,21),...
                   0.2*(omega(4)-omega(3))*linspace(-1,1,21));
  save(filename,'w1','w2','DM');  
else
  clear DM*
  load(filename)
end;

for k=1:length(DM),
  variable = ['DM',DM{k}];
  disp([variable,'=dm;']);

  var = whos('-file',filename);
  j = find(strcmp({var(:).name},variable)==1);

  if isempty(j),
    fprintf('============== %s ====================\n\n',variable)
    dm = zeros(size(w1));
    for j=1:numel(w1),
      Y = trafo([w1(j);w2(j)],xc(:));
      Tc = inter(T,omega,Y);
      dm(j) = feval(DM{k},Tc,Rin,omega,m,'edge',edge);
      if j== 1,
        figure(3);
        ph = viewImage(Tc,omega,m);
      else
        set(ph,'cdata',reshape(Tc,m)'); drawnow
        pause(1/100)
      end;
      title(sprintf('%d/%d:%d/%d',k,length(DM),j,numel(w1)));
    end;
    eval([variable,'=dm;']);
    save(filename,'-append',variable);
  end;
end;

load(filename)

%figure(1); close(1); figure(1); clf; set(1,'position',position(800),'color','w');
%figure(2); close(2); figure(2); clf; set(2,'position',position(800),'color','w');
meshPrint = @(k) FAIRprint(sprintf('%s-mesh',DM{k}));
contPrint = @(k) FAIRprint(sprintf('%s-contour',DM{k}));
shift = [3000,0.5,-1,0.08];
Ztick = {[0:1000:4000],[0:0.1:0.3],[-1:0:1],[0:0.02:0.06]};

for k=1:4,

  variable = ['DM',DM{k}];
  eval(['dm=',variable,';']);

  fig = figure(10+k);
  close(fig); figure(fig); set(fig,'position',FAIRposition(800),'color','w');
  ph=mesh(w1,w2,dm-shift(k)); grid off;
  set(gca,'fontsize',30);
  view(-135,25);
  set(ph,'linewidth',2)
  meshPrint(k);
  
  fig = figure(20+k);
  close(fig); figure(fig); clf; set(fig,'position',FAIRposition(800),'color','w');
  contour(w1,w2,dm,10,'linewidth',2);
  set(gca,'fontsize',30);
  contPrint(k);
end;
