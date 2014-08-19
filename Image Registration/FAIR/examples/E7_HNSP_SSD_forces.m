%==============================================================================
% (c) Jan Modersitzki 2011/01/02, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% Tutorial for FAIR: SSD and Force Fields
%
%==============================================================================

clear, close all, help(mfilename);
sdiag = @(a) spdiags(reshape(a,[],1),0,length(a),length(a));
setupHNSPData; 
level = 7; omega = MLdata{level}.omega; m = MLdata{level}.m; % load data
inter('reset','inter','splineInter')
viewImage('set','axis','off');
[T,R] = inter('coefficients',MLdata{level}.T,MLdata{level}.R,omega);
X = reshape(getCellCenteredGrid(omega,m),[],2);

[Tc,dT] = inter(T,omega,X);
[Rc,dR] = inter(R,omega,X);

axT   = [0.55,0.85,0.32,0.47];
axR   = [0.6,0.9,0.20,0.35];
axTF  = [0.55,0.70,0.32,0.395];

FAIRprint('reset','folder',fullfile(FAIRpath,'temp'),'prefix','HNSP-',...
  'pause',1,'obj','gca','format','jpg','draft','on');
FAIRprint('disp');

FAIRfigure(1,'color','w','position',800);
viewImage(Tc,omega,m); hold on; axis off;
ph = plot(axT([1,1,2,2,1]),axT([3,4,4,3,3]),'w-','linewidth',4); 
FAIRprint('T0');

set(ph,'visible','off');  axis(axT);  FAIRprint('T0-ROI');

clf; 
viewImage(Rc,omega,m); hold on; axis off;
ph = plot(axR([1,1,2,2,1]),axR([3,4,4,3,3]),'w-','linewidth',4);
FAIRprint('R0')

set(ph,'visible','off');  axis(axR);  FAIRprint('R0-ROI');

clf;
viewImage(128+0.5*(Tc-Rc),omega,m); hold on; axis off;
% ph = plot(axT([1,1,2,2,1]),axT([3,4,4,3,3]),'w-','linewidth',4);
FAIRprint('D0')
% set(ph,'visible','off');  axis(axR);  FAIRprint('D0-ROI')

X = getCellCenteredGrid(omega,m); X = reshape(X,[],2);
n = prod(m);
gradT     = spdiags(dT,[0,n]);      
maxGradT  = norm(gradT,'inf');
ngradT    = 128+128/maxGradT*gradT;
normGradT = sqrt(sum(gradT.^2,2));
J     = find(normGradT>5e1);
F     = sdiag(Tc-Rc)*gradT;      maxF = norm(F,'inf');
nF    = 128+128/maxF*F;
normF = sqrt(sum(F.^2,2));
K     = find(normF>1e2);

figure(1); clf; set(1,'position',FAIRposition(800),'color','w');
clf; viewImage(ngradT(:,1),omega,m); axis(axT); FAIRprint('d1T');
clf; viewImage(ngradT(:,2),omega,m); axis(axT); FAIRprint('d2T');

clf; viewImage(128+0.5*Tc,omega,m); hold on;
qh = quiver(X(J,1),X(J,2),gradT(J,1),gradT(J,2),2);
set(qh,'linewidth',3,'color','k'); axis(axTF);
FAIRprint('gradT');

clf; viewImage(128+0.5*Tc,omega,m); hold on;
qh = quiver(X(K,1),X(K,2),F(K,1),F(K,2),2);
set(qh,'linewidth',3,'color','k'); axis(axTF)
FAIRprint('F');
