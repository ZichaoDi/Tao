%% Analytical formula to find the shift of each projection due to center of rotation transformation;
delta_d=0; % off center for the initial reference projection;
Det=norm(DetKnot0(1,:)-SourceKnot0(1,:));
Num=(DetKnot0(1,1)-SourceKnot0(1,1))*((cos(theta')-1).*delta0(:,2)'+delta0(:,1)'.*sin(theta'))+ (DetKnot0(1,2)-SourceKnot0(1,2))*((1-cos(theta')).*delta0(:,1)'+sin(theta').*delta0(:,2)');
% Mt=XRF_decom(:,:,5)';
shift=shift(1,:);%xCOR(1:N_delta);
% shift=Num./Det/dTau+1/2*delta_d0/dTau;
% Mt=-log(DisR./I0');
alignedDiscrete=zeros(nTau+1,numThetan);
alignedContinuous=alignedDiscrete;
Mt=b;
tt=Mt;%Mt(:,:,2);
for i = 1:numThetan
    delay=shift(i);
    %%---------------- Continuous: Gaussian Convolution
    sigma=1.5/2.355;
    nT=nTau+1;
    %%=============================================
    H=fft(tt(:,i)');
    u=[0:floor((nT)/2)-1 floor(-(nT)/2):-1] / (nT);
    ubar=exp(-sqrt(-1)*2*pi.*(u*delay));
    H1=H.*ubar;
    alignedContinuous(:,i)=real(ifft(H1));
    %%=============================================
    H=fftshift(fft(tt(:,i)'));
    u=[1:nTau+1]/(nTau+1);
    ubar=exp(-sqrt(-1)*2*pi.*(u*delay));
    H1=H.*ubar;
    alignedContinuous1(:,i)=(abs(ifft((H1))));
    %%=========================================================
    scale=1/(sqrt(2*pi)*sigma);
    range=[0:floor((nTau+1)/2)-1 floor(-(nTau+1)/2):-1]';% / (nTau+1);
    G=exp(-(range-delay).^2./(2*sigma^2));
    alignedContinuous2(:,i)=scale*real(ifft((fft(G)).*(fft(tt(:,i)))));
    %%---------------- Discrete: Permutation Matrix
    if(delay>=0)
        alignedDiscrete(delay+1:end,i)=tt(1:end-delay,i);
        alignedDiscrete(1:delay,i)=tt(end-delay+1:end,i);
    else
        alignedDiscrete(1:end+delay,i)=tt(-delay+1:end,i);
        alignedDiscrete(end+delay+1:end,i)=tt(1:-delay,i);
    end
    %%--------------------------------------------------
    subplot(1,2,1)
    plot(Mt(:,i,1),'r.-');
    hold on;
    plot(alignedContinuous(:,i),'m.--');
    plot(alignedContinuous2(:,i),'k.-');
    title(num2str([delay,i]));
    hold off;
    subplot(1,2,2)
    plot(range,G,'r.-');
    hold off;
    pause;
end
figure, 
nrow=4;
subplot(nrow,2,1);imagesc(tt); title('Observed Sinogram')
subplot(nrow,2,3);imagesc(alignedContinuous2);title('Continuous Corrected')
% subplot(nrow,2,4);imagesc(iradon(fliplr(alignedContinuous),thetan)');%,'linear','shepp-logan',N)'); 
% axis xy image;
subplot(nrow,2,5);imagesc(alignedDiscrete);title('Discrete Corrected')
% subplot(nrow,2,6);imagesc(iradon(fliplr(alignedDiscrete),thetan)');%,'linear','shepp-logan',N)'); 
% axis xy image;
subplot(nrow,2,7);imagesc(Mt(:,:,1));title('Sinogram with COR(0,0)')
subplot(nrow,2,[2 4 6 8]);imagesc(W);axis xy image;hold on;
plot(delta0(:,1)/dz(1)+N/2,delta0(:,2)/dz(2)+N/2,'r.-');
