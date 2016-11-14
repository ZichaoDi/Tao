load PeriodicTable
% Z=[79 14 74];
% NumElement=3;
% load tomopytest.mat
% N=195;
% a=data_xrf_decom;
% figure, 
% for i=1:3; 
%     subplot(1,3,i);imagesc(map1D(squeeze(a(i,:,:)),[0 1]).^0.4);set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]); axis xy image; 
%     title(Element{Z(i)},prop_names,prop_values);
%     if(i==1);
%         xlabel('$\tau$',prop_names,prop_values);
%         ylabel('$\theta$',prop_names,prop_values);
%     end
% end
% h=colorbar('SouthOutside'); 
% set(h, 'Position', [0.45 .4 .4 .02]);
% print(figtype,'../figures/decom_xrf')
LW = 2; % Linewidth
FS = 14; % Fontsize
MS = 10; % MarkerSize

marker={'o','v','s','d'};
figtype = '-depsc';

prop_names(1) = {'FontWeight'};
prop_names(2) = {'FontSize'};
prop_names(3) = {'Interpreter'};
prop_values(1,1) = {'bold'};
prop_values(1,2) = {FS};
prop_values(1,3) = {'latex'};

load spectra_rod_downsample
% load spectra_30_aligned
% load backgroud_subtract
% BG=reshape(spectra_30_aligned,1748*73,2000);
% BG=sum(BG(ind_bg,:),1)'*1;
load BerylliumTransmission
load simulate_recon195;
b=squeeze(sum(sum(XRF,1),2)).*DetMu;
figure,
h1=semilogy(DetChannel,squeeze(a),'k-','LineWidth',2);
hold on;
h2=semilogy(DetChannel,b/0.19+1,'m-','LineWidth',2);%+BG*4.2  
NumElement=8;
cmap=colormap(colorcube(NumElement+1));%rand(NumElement,3);
cmap=cmap([1:NumElement],:);
p=[];
ub=max(a(:))+1e4;
maxE=13;
for ele=1:NumElement
    for el =1:size(BindingEnergy,2)
        if(BindingEnergy(ele,el)>0 & BindingEnergy(ele,el)<maxE)

             p{ele}=plot([BindingEnergy(ele,el) BindingEnergy(ele,el)],[1 ub],'Color',cmap(ele,:),'LineStyle','--');
             hold on;
             % text(double(BindingEnergy(ele,el)),double(ub),Element(Z(ele)));
          end
   end
end
xlabel('Energy (keV)',prop_names,prop_values);
ylabel('Intensity (counts/sec)',prop_names,prop_values);
% l=legend([p{1},p{2},p{3},p{4},p{5},p{6},p{7},p{8},h1,h2],{'Au','B','O','Na','Al','Si','K','W','experimental spectrum','simulated spectrum on reconstructed sample'},'Location','best');
% % l=legend([p{1},p{2},p{3},h1,h2],{'Au','Si','W','experimental spectrum','simulated spectrum on reconstructed map'},'Location','best');
% title(l,'Spectra and Emission Lines');
axis([0 maxE 0 ub])
% print(figtype,'../figures/spec_rod_sum')
