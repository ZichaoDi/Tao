%% Analytical formula to find the shift of each projection due to center of rotation transformation;
delta_d=0; % off center for the initial reference projection;
Det=norm(DetKnot0(1,:)-SourceKnot0(1,:));
Num=(DetKnot0(1,1)-SourceKnot0(1,1))*((cos(theta')-1).*delta0(:,2)'+delta0(:,1)'.*sin(theta'))+ (DetKnot0(1,2)-SourceKnot0(1,2))*((1-cos(theta')).*delta0(:,1)'+sin(theta').*delta0(:,2)');
shift=Num./Det/dTau+1/2*delta_d0/dTau;
Mt=-log(DisR./I0');
alignedDiscrete=zeros(nTau+1,numThetan);
alignedContinuous=alignedDiscrete;
alignedContinuous1=alignedDiscrete;
alignedContinuous2=alignedDiscrete;
Mt1=Mt(:,:,1);%Mth(1:4/factor:end,:,1);
tt=Mt(:,:,2);
sig=linspace(0.0125,1.0616,50);
err=[];
shift=15*ones(numThetan,1);
for ti=1:length(sig)
    sigma=sig(ti);
for i = 1:numThetan
    delay=shift(i);
    %%---------------- Continuous: Gaussian Convolution
    % sigma=1.5/2.355;
    nT=nTau+1;
    %%=============================================
    %%=============================================
    % H=fftshift(fft(tt(i,:)));
    % u=[1:nTau+1]/(nTau+1);
    % ubar=exp(-sqrt(-1)*2*pi.*(u*delay));
    % H1=H.*ubar;
    % alignedContinuous1(:,i)=(abs(ifft((H1))));
    %%=========================================================
    scale=1/(sqrt(2*pi)*sigma);
    range=[0:floor((nTau+1)/2)-1 floor(-(nTau+1)/2):-1]';% / (nTau+1);
    G=exp(-(range-delay).^2./(2*sigma^2));
    alignedContinuous(:,i)=scale*real(ifft((fft(G)).*(fft(tt(i,:)'))));
    %%---------------- Discrete: Permutation Matrix
    delay=round(delay);
    if(delay>=0)
        alignedDiscrete(delay+1:end,i)=tt(i,1:end-delay);
        alignedDiscrete(1:delay,i)=tt(i,end-delay+1:end);
    else
        alignedDiscrete(1:end+delay,i)=tt(i,-delay+1:end);
        alignedDiscrete(end+delay+1:end,i)=tt(i,1:-delay);
    end
    %%--------------------------------------------------
    % subplot(1,2,1)
    % plot(Mt(:,i,1),'r.-');
    % hold on;
    % plot(alignedContinuous(:,i),'m.--');
    % plot(alignedContinuous2(:,i),'k.-');
    % title(num2str([delay,i]));
    % hold off;
    % subplot(1,2,2)
    % plot(range,G,'r.-');
    % hold off;
    % pause;
end
% err(ti,:)=[norm(Mt1(:)-alignedContinuous(:)),norm(Mt1(:)-alignedDiscrete(:))];
err(ti,:,:)=Mt1-alignedContinuous';
end
% save('/nfs2/wendydi/Documents/Research/optdi/NanoLDRD/SISC/GenerateFig/data/CoRanalytic_1mri.mat','Mt1','Mt2','delta0','dz','N','alignedDiscrete','alignedContinuous','W')
figure, 
nrow=3;
subplot(nrow,2,1);imagesc(Mt(:,:,2)); title('Observed Sinogram')
xlabel('\theta');ylabel('\tau','FontWeight','bold','FontSize',14);
% set(gca,'xtick',[],'ytick',[]);
% subplot(nrow,2,2);imagesc(reshape(W,N,N));%,'linear','shepp-logan',N)'); 
subplot(nrow,2,2);imagesc(iradon(Mt(:,:,2),thetan));%,'linear','shepp-logan',N)'); 
axis xy image;
set(gca,'xtick',[],'ytick',[]);
title('Reconstruction w/o correction')
subplot(nrow,2,3);imagesc(alignedContinuous');title('Continuous Corrected')
% set(gca,'xtick',[],'ytick',[]);
subplot(nrow,2,4);imagesc(reshape(W,N,N),[0,max(W(:))]);%,'linear','shepp-logan',N)'); 
title('Reconstruction with correction')
set(gca,'xtick',[],'ytick',[]);
axis xy image;
% subplot(nrow,2,5);imagesc(alignedDiscrete);title('Discrete Corrected')
% subplot(nrow,2,6);imagesc(iradon(fliplr(Mt(:,:,2)),thetan)');%,'linear','shepp-logan',N)'); 
% axis xy image;
subplot(nrow,2,5);imagesc(Mt(:,:,1));title('Sinogram with COR(0,0)')
% set(gca,'xtick',[],'ytick',[]);
subplot(nrow,2,6);
imagesc(W);axis xy image;hold on;
set(gca,'xtick',[],'ytick',[]);
title('Ground Truth')
plot(delta0(:,1)/dz(1)+N/2,delta0(:,2)/dz(2)+N/2,'r*-','MarkerSize',5);
return;
subplot(1,2,1);imagesc(Mt(:,:,1));
set(gca,'xtick',[],'ytick',[]);
xlabel('\theta');ylabel('\tau','FontWeight','bold','FontSize',14);
title('Sinogram w/o CoR shifts')
subplot(1,2,2);imagesc(Mt(:,:,2));
set(gca,'xtick',[],'ytick',[]);
xlabel('\theta');ylabel('\tau','FontWeight','bold','FontSize',14);
title('Sinogram with CoR shifts')
