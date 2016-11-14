close all
LW = 2; % Linewidth
FS = 18; % Fontsize
MS = 10; % MarkerSize

marker={'o','v','s','d'};
figtype = '-depsc';
figPath = '../figures/';
prop_names(1) = {'FontWeight'};
prop_names(2) = {'FontSize'};
prop_names(3) = {'Interpreter'};
prop_values(1,1) = {'bold'};
prop_values(1,2) = {FS-2};
prop_values(1,3) = {'Tex'};

load /nfs2/wendydi/Documents/Research/Tao/XRF_XTM_Simulation/XRF_Simpler/data_rod
N=160;
load PeriodicTable
Z=[79 14 74];
NumElement=length(Z);
xsR=reshape(xstar_xrf,N,N,NumElement);
xsJ=reshape(xstar_joint,N,N,NumElement);;
xsT=reshape(xstar_tomopy,195,195,NumElement);;
w0=reshape(x0,N,N,NumElement);;
numRow=4;
leftF=0.1;
bottomF=0.7;
widthF=0.15;
lengthF=0.15;
positionVector1 = [leftF bottomF widthF lengthF];
shift_l=[1.1 0 0 0];
shift_b=[0 1.1 0 0];
for i=1:NumElement,
    subplot('position',positionVector1+shift_l*widthF*(i-1)),
    imagesc(w0(:,:,i));
    if(i==1)
        ylabel({'Initial';'Guess'},prop_names,prop_values)
    end
    title(Element{Z(i)},prop_names,prop_values);
    colormap(jet),
    set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]);
    
end

positionVector1 = positionVector1-shift_b*lengthF;
for i=1:NumElement,
    subplot('position',positionVector1+shift_l*widthF*(i-1)),
    imagesc(xsR(:,:,i));
    if(i==1)
        ylabel('XRF',prop_names,prop_values)
    end
    colormap(jet),
    set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]);
end
positionVector1 = positionVector1-shift_b*lengthF;
for i=1:NumElement,
    subplot('position',positionVector1+shift_l*widthF*(i-1)),
    imagesc(xsJ(:,:,i));
    if(i==1)
        ylabel('Joint',prop_names,prop_values)
    end
    colormap(jet),
    set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]);
end
h=colorbar('eastOutside');
 positionVector1 = positionVector1-shift_b*lengthF;
 for i=1:NumElement,
     subplot('position',positionVector1+shift_l*widthF*(i-1)),
     imagesc(xsT(:,:,i)');
     if(i==1)
         ylabel('TomoPy',prop_names,prop_values)
     end
     colormap(jet),
     set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]);
 end
% set(h, 'Position', [.1 .05 .30 .02]);%
set(h, 'Position', [0.6 .41 .01 .4]);
% print(figtype,'../figures/RodResult')
