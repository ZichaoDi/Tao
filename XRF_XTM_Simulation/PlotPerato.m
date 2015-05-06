FS = 16;
MS = 10;
LW = 1.5;
marker={'o','v','s','d'};

load ./result/MatFile/betaPerf_phantom3_4.mat
%    load ./result/MatFile/betaPerf_phase20_12.mat
[m(1),n(1)]=size(Abeta2);
[m(2),n(2)]=size(Abeta5);
[m(3),n(3)]=size(AbetaInf);

Iter=[1,16]; %% Iteration numbers particular to compare
hl=[];LegendInf=[];
close all
hl(1)=loglog(Abeta5(:,5),Abeta5(:,6),'b-');
LegendInf{1}='JRT';
hold on;
hl(2)=loglog(Abeta2(:,5),Abeta2(:,6),'r-.');
LegendInf{2}='XRF';
hl(3)=loglog(AbetaInf(:,5),AbetaInf(:,6),'m--');
LegendInf{3}='XRT';

 %% plot error reduction with respect to XRF and XRT
for i=1:length(Iter)
    hl(i+3) = loglog(Abeta5(Iter(i),5),Abeta5(Iter(i),6),['k',marker{i}]);
    loglog(Abeta5(Iter(i),5),Abeta5(Iter(i),6),['b',marker{i}],'MarkerFaceColor','b');
    loglog(Abeta2(Iter(i),5),Abeta2(Iter(i),6),['r',marker{i}],'MarkerFaceColor','r');
    loglog(AbetaInf(Iter(i),5),AbetaInf(Iter(i),6),['m',marker{i}],'MarkerFaceColor','m');
    LegendInf{i+3}=['Iteration# ',num2str(Iter(i)-1)];
end
%%---------If plot the end iteration
% hl(6)=loglog(Abeta5(end,5),Abeta5(end,6),['b',marker{3}],'MarkerFaceColor','b');
% LegendInf{6}=['Iteration# ',num2str(m(2)-1)];
% hl(7)=loglog(Abeta2(end,5),Abeta2(end,6),['r',marker{4}],'MarkerFaceColor','r');
% %     loglog(AbetaInf(Iter(i),5),AbetaInf(Iter(i),6),['m',marker{i}],'MarkerFaceColor','m');
% LegendInf{7}=['Iteration# ',num2str(m(1)-1)];
%%-----------------------------------------------------
hleg=legend(hl,LegendInf);
set(findall(gcf,'type','line'),'LineWidth',LW,'MarkerSize',MS)
set(hleg,'FontSize',FS,'FontWeight','bold','Location','Best');
xlabel('XRT Residual (Log scaled)','FontSize',FS,'FontWeight','bold')
ylabel('XRF Residual (Log-scaled)','FontSize',FS,'FontWeight','bold')
set(gca,'FontSize',FS-2)
axis tight
box on
% print('-depsc', ['~/Documents/TAO/TAOlocal/optdi/Reconstruction/Writing/figures/Residual_p',num2str(m(1)),'_',num2str(numTheta)]);
%%=================================================================================
figure, %% Plot residual reduction v.s. f-g evaluations
hr=semilogy(1:size(Abeta5,1),Abeta5(:,end),'bv-',1:size(Abeta2,1),Abeta2(:,end),'rs-');
hleg=legend(hr,LegendInf{1:2});
set(findall(gcf,'type','line'),'LineWidth',LW,'MarkerSize',MS)
set(hleg,'FontSize',FS,'FontWeight','bold','Location','Best');%'southeast');
xlabel('Number of f-g Evaluations','FontSize',FS,'FontWeight','bold')
ylabel('Reconstruction Error (Log-scaled)','FontSize',FS,'FontWeight','bold')
set(gca,'FontSize',FS-2)
axis tight
box on
% print('-depsc', '~/Documents/TAO/TAOlocal/optdi/Reconstruction/Writing/figures/Error_p3_4')
