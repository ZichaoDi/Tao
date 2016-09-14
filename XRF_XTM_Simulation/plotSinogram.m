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

% load spectra_rod
% load resulted_spectra1
load simulate_realrod_solid
b=map1D(squeeze(sum(sum(XRF,1),2)),a)+1;
figure,
h1=semilogy(DetChannel,squeeze(a),'k-','LineWidth',1);
hold on;
h2=semilogy(DetChannel,b,'m-','LineWidth',1);  
cmap=rand(NumElement,3);%[1 0 0; 0 1 0; 0 0 1];
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
xlabel('Energy (keV)');%,prop_names,prop_values);
ylabel('Intensity (counts/sec)');%,prop_names,prop_values);
l=legend([p{1},p{2},p{3},p{4},p{5},p{6},p{7},h1,h2],{'Au','B','O','Na','Al','Si','W','experimental spectrum','simulated spectrum on real data'},'Location','best');
title(l,'Spectra and Emission Lines');
axis([0 maxE 0 ub])
% print(figtype,'../figures/spec_rod_sum')
