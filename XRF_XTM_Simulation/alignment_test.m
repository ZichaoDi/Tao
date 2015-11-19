% P=phantom(128);
P=checkerboard(12);
theta=linspace(0,180,180);
R=radon(P,theta);
rng('default')
R_shift=R;
rand_shift=floor(0+rand(size(R,2),1)*8);
for i=2:2:size(R,2);
    R_shift(:,i)=circshift(R(:,i),rand_shift(i));
end
probe_drift=0+10*rand(size(theta));
theta_shift=theta+probe_drift;
R_drift=radon(P,theta_shift);
xt=linspace(0,1,size(R,1))';
drift=1/(sigma*sqrt(2*pi))*exp(-(xt-miu).^2./(2*sigma^2));%sort(18*rand(size(R(:,1))));
R_drift=R_drift+repmat(10*drift,1,size(R_drift,2));
R_s=align_center(R_shift);
P_shift=iradon(R_shift,theta);
P_s=iradon(R_s,theta);
P_d=iradon(R_drift,theta);
P0=iradon(R,theta);
colNum=2;
figure(10), subplot(2,colNum,1),imagesc(R),title('Perfect Sinogram')
set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]);
% subplot(2,colNum,2),imagesc(R_shift),title('Misaligned Sinogram')
% set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]);
subplot(2,colNum,2),imagesc(R_drift),title('probe drift Sinogram')
set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]);
%subplot(2,colNum,3),imagesc(R_s);title('Aligned Sinogram')
%set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]);
subplot(2,colNum,3),imagesc(P);title('True Sample')
%set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]);
%subplot(2,colNum,6),imagesc(P_shift);title({'Reconstruction';' without alignment'})
%set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]);
%subplot(2,colNum,7),imagesc(P_s);title({'Reconstruction';'with alignment'})
%set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]);
subplot(2,colNum,4),imagesc(P_d);title({'Reconstruction probe drift'})
set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]);

    
    
