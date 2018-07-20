%% Analytical formula to find the shift of each projection due to center of rotation transformation;
delta_d=0; % off center for the initial reference projection;
Mt=-log(DisR./I0');
alignedDiscrete=zeros(nTau+1,numThetan);
alignedContinuous=alignedDiscrete;
tt=Mt(:,:,2);
shift=pert_scan*nTau;
for i = 1:numThetan
    delay=shift(i);
    %%---------------- Continuous: Gaussian Convolution
    sigma=1.5/2.355;
    nT=nTau+1;
    %%=============================================
    scale=1/(sqrt(2*pi)*sigma);
    range=[0:floor((nTau+1)/2)-1 floor(-(nTau+1)/2):-1]';% / (nTau+1);
    G=exp(-(range-shift).^2./(2*sigma^2));
    alignedContinuous(:,i)=scale*real(ifft((fft(G)).*(fft(tt(:,i)))));
    pause;
    %%---------------- Discrete: Permutation Matrix
    if(delay>=0)
        alignedDiscrete(delay+1:end,i)=tt(1:end-delay,i);
        alignedDiscrete(1:delay,i)=tt(end-delay+1:end,i);
    else
        alignedDiscrete(1:end+delay,i)=tt(-delay+1:end,i);
        alignedDiscrete(end+delay+1:end,i)=tt(1:-delay,i);
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
figure, 
nrow=3;
subplot(nrow,2,1);imagesc(Mt(:,:,2)); title('Observed Sinogram')
xlabel('\theta');ylabel('\tau','FontWeight','bold','FontSize',14);
% set(gca,'xtick',[],'ytick',[]);
subplot(nrow,2,2);imagesc(reshape(W,N,N));%,'linear','shepp-logan',N)'); 
axis xy image;
set(gca,'xtick',[],'ytick',[]);
title('Reconstruction w/o correction')
subplot(nrow,2,3);imagesc(alignedContinuous);title('Continuous Corrected')
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
return;
plot(delta0(:,1)/dz(1)+N/2,delta0(:,2)/dz(2)+N/2,'r*-','MarkerSize',5);
subplot(1,2,1);imagesc(Mt(:,:,1));
set(gca,'xtick',[],'ytick',[]);
xlabel('\theta');ylabel('\tau','FontWeight','bold','FontSize',14);
title('Sinogram w/o CoR shifts')
subplot(1,2,2);imagesc(Mt(:,:,2));
set(gca,'xtick',[],'ytick',[]);
xlabel('\theta');ylabel('\tau','FontWeight','bold','FontSize',14);
title('Sinogram with CoR shifts')
