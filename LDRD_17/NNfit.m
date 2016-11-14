Tomo_startup;
load PeriodicTable
load UnitSpectrum
load DetChannel_Rod
element={'Si', 'Cl', 'Ar', 'K', 'Ca', 'Ti', 'Fe', 'Cu', 'Zn', 'W_L', 'Au_L','Au_M', 'inelastic_scatter', 'elastic_scatter'};
Z=[14 17 18 19 20 22 26 29 30 74 79];
nrow=2000;
ncol=14;
M_sv=zeros(4,nrow,ncol);
for i=1:4;
    M_sv(i,:,:)=double(csvread(['specs_det_',num2str(i-1),'.csv'],0,0,[0 0 nrow-1 ncol-1]));
end

load spectra_30_aligned;
d=double(squeeze(sum(sum(spectra_30_aligned,1),2)));
ub=max(d);
x=zeros(4,ncol);
for i=1:4
    x(i,:) = lsqnonneg(squeeze(M_sv(i,:,:)),d);
end
x_unit=lsqnonneg(M_raw(Z,:)',d);
fitted_d_unit= M_raw(Z,:)'*x_unit;
%%================plot fitted spectrum
for ind=1:1;
    fitted_d= squeeze(M_sv(ind,:,:))*x(ind,:)';
    figure,subplot(2,1,1);semilogy(DetChannel,d,'k.-',DetChannel,fitted_d,'m.--');
    axis([0 13 1 ub]);
    legend('Experimental Spectrum','Fitted Spectrum','Location','best')
    title(['specs-det',num2str(ind)]);
    subplot(2,1,2);plot(1:ncol, x(ind,:),'r.--');
    xlabel('index of known spectra')
    ylabel({'weight of';'each known spectrum'})
    for i=1:ncol, text(i,x(ind,i),element{i});end
end

figure,subplot(2,1,1);semilogy(DetChannel,d,'k.-',DetChannel,fitted_d_unit,'m.--');
axis([0 13 1 ub]);
legend('Experimental Spectrum','Fitted Spectrum','Location','best')
subplot(2,1,2);plot(1:length(Z), x_unit,'r.--');
for i=1:length(Z), text(i,x_unit(i),Element(Z(i)));end

           
