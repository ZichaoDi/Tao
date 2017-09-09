%% Analytical formula to find the shift of each projection due to center of rotation transformation;
% c=15;
% delta0(1)= c*Tau/1748;
% delta0(2)= c*Tau/1748;
global Ddelta0_COR

delta_d=0; % off center for the initial reference projection;
step=[0];%-10:1e-1:10;
eps1=zeros(numThetan,1);

Det=norm(DetKnot0(1,:)-SourceKnot0(1,:));

eps2_Det=Det*dTau;
eps2_Num=...
[2*(DetKnot0(1,1)-SourceKnot0(1,1))*(-sin(theta')+sin(2*theta'))+2*(DetKnot0(1,2)-SourceKnot0(1,2))*(cos(2*theta')-cos(theta')); ...
2*(DetKnot0(1,1)-SourceKnot0(1,1))*(cos(theta')-cos(2*theta'))+2*(DetKnot0(1,2)-SourceKnot0(1,2))*(sin(2*theta')-sin(theta'))];
eps2=eps2_Det./eps2_Num./dz(1);
eps2(isinf(eps2))=0;

Ddelta0=[(DetKnot0(1,1)-SourceKnot0(1,1))*(-sin(theta')+sin(2*theta'))+(DetKnot0(1,2)-SourceKnot0(1,2))*(cos(2*theta')-cos(theta'));...
(DetKnot0(1,1)-SourceKnot0(1,1))*(cos(theta')-cos(2*theta'))+(DetKnot0(1,2)-SourceKnot0(1,2))*(sin(2*theta')-sin(theta'))];
Ddelta0=Ddelta0/eps2_Det;
Ddelta0_COR=repmat(Ddelta0,[1 nTau+1]);
k_outer=[1];
err_step=zeros(length(k_outer),length(step));
shift=zeros(length(step),numThetan);
for ko_ind=1:length(k_outer);

for k=1:length(step)
    k_step=step(k); 
    Num=(DetKnot0(1,1)-SourceKnot0(1,1))*((cos(theta')-1).*delta0(:,2)'+delta0(:,1)'.*sin(theta'))+ (DetKnot0(1,2)-SourceKnot0(1,2))*((1-cos(theta')).*delta0(:,1)'+sin(theta').*delta0(:,2)');
    % Mt=XRF_decom(:,:,2)';
    shift(k,:)=Num./Det/dTau+1/2*delta_d0/dTau;
    Mt=-log(DisR./I0');
    alignedDiscrete=zeros(nTau+1,numThetan);
    LR2=zeros(nTau+1,numThetan);
    alignedContinuous=alignedDiscrete;
    tt=Mt(:,:,2);
    nfill=0;;
    for i = 1:numThetan
        % delay=floor(mod(shift(k,i),nTau+1));
        % delay=double(mod(int32(shift(k,i)),nTau+1));
        delay=shift(k,i);
        %%---------------- Continuous: Gaussian Convolution
        sigma=0.001/2.355;
        nT=nTau+1;
        %%=============================================
        % nT=1000;
        % nTau=nT;
        % G1=exp(-([1:nTau+1]'-80).^2./(2*1^2));
        % G2=exp(-([1:nTau+1]'-200).^2./(2*2^2));
        % G3=exp(-([1:nTau+1]'-500).^2./(2*10.6^2));
        % G4=exp(-([1:nTau+1]'-700).^2./(2*0.5^2));
        % Gt=zeros(nTau+1,1);Gt(25:31)=1;
        % fine_grid=1:1:nT;
        % tzip=G1+G2+G3+G4+Gt;
        % nT=length(tzip);
        % delay=0.24;
        % H=fft(tzip');
        % u=[0:floor((nT)/2)-1 floor(-(nT)/2):-1] / (nT);
        % ubar=exp(-sqrt(-1)*2*pi.*(u*delay));
        % H1=H.*ubar;
        % alignedContinuous=(real(ifft(H1)));
        % figure,plot(tzip,'r.-');hold on; plot(alignedContinuous,'b.-')
        % return;
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
        % uu=interp1([1:nTau+1],tt(:,i),[1:1:nTau+1]);
        % H=fft(tt(:,i));
        % u=[0:floor((nTau+1)/2)-1 floor(-(nTau+1)/2):-1]' / (nTau+1);
        % u=exp(-sqrt(-1)*2*pi.*(u*delay));
        % H1=H.*u;
        % alignedContinuous1(:,i)=real(ifft(H1));
        %%=========================================================
        G=exp(-([-nTau-1:nTau+1]'-delay).^2./(2*sigma^2));
        aligned_temp=2/pi*real(ifft(fft(G).*fft([0*ones(nTau+2,1);tt(:,i)])));
        if(abs(delay)>(nTau+1)/2)
            alignedContinuous3(:,i)=aligned_temp(nTau+2:end-1);
        else
            alignedContinuous3(:,i)=aligned_temp(1:nTau+1);
        end
        %%=========================================================
        range=[0:floor((nTau+1)/2)-1 floor(-(nTau+1)/2):-1]';% / (nTau+1);
        G=exp(-(range-delay).^2./(2*sigma^2));
        alignedContinuous2(:,i)=2/2*real(ifft((fft(G)).*(fft(tt(:,i)))));
        %%---------------- Discrete: Permutation Matrix
        % if(delay>=0)
        %     alignedDiscrete(delay+1:end,i)=tt(1:end-delay,i);
        %     alignedDiscrete(1:delay,i)=tt(end-delay+1:end,i);
        % else
        %     alignedDiscrete(1:end+delay,i)=tt(-delay+1:end,i);
        %     alignedDiscrete(end+delay+1:end,i)=tt(1:-delay,i);
        % end
        %%--------------------------------------------------
        % plot(Mt(:,i,1),'r.-');
        % hold on;
        % plot(Mt(:,i,2),'m.--');
        % % plot(alignedContinuous(:,i),'b.-');title(num2str([delay,i]));
        % plot(alignedContinuous2(:,i),'k.-');title(num2str([delay,i]));
        % % plot(alignedContinuous1(:,i),'g.-');title(num2str([delay,i]));
        % plot(alignedContinuous3(:,i),'b.-');title(num2str([delay,i]));hold off;
        % hold off;
        % pause;
        % subplot(1,3,2),plot(Mt(:,i,1),'r.-');hold on; plot(alignedContinuous1(:,i),'b.-');title(num2str(delay));hold off;
        % subplot(1,3,3),plot(Mt(:,i,1),'r.-');hold on; plot(alignedDiscrete(:,i),'b.-');title(num2str(delay));hold off;
        % pause;
    end
    %%---------------- Continuous: Gaussian Convolution
    % G=1/(sigma*sqrt(2*pi))*exp(-(repmat([1:nTau+1]',1,numThetan)-repmat(mod(shift(k,:),nTau+1),nTau+1,1)).^2./(2*sigma^2));
    % alignedContinuous1=ifft(fft(G).*fft(tt));
    % errD=alignedContinuous1-alignedContinuous;
    %%--------------------------------------------------
end
end
figure, 
nrow=4;
subplot(nrow,2,1);imagesc(tt); title('Observed Sinogram')
% axis xy image;
% subplot(nrow,2,2);imagesc(iradon(fliplr(tt),thetan,'linear','shepp-logan',N)'); 
% axis xy image;
subplot(nrow,2,3);imagesc(alignedContinuous);title('Continuous Corrected')
% subplot(nrow,2,4);imagesc(iradon(fliplr(alignedContinuous),thetan)');%,'linear','shepp-logan',N)'); 
% axis xy image;
% subplot(nrow,2,5);imagesc(alignedDiscrete);title('Discrete Corrected')
% subplot(nrow,2,6);imagesc(iradon(fliplr(alignedDiscrete),thetan)');%,'linear','shepp-logan',N)'); 
% axis xy image;
subplot(nrow,2,7);imagesc(Mt(:,:,1));title('Sinogram with COR(0,0)')
% axis xy image;
subplot(nrow,2,[2 4 6 8]);imagesc(W);axis xy image;hold on;
% plot(delta0(:,1)/dz(1)+N/2,delta0(:,2)/dz(2)+N/2,'r.-'); 
% for i=1:numThetan
%     text(delta0(i,1)/dz(1)+N/2,delta0(i,2)/dz(2)+N/2,num2str(i));
% end
% axis xy image;


% d1=[];d2=[];for i=1:numThetan;[~,d]=max(Mt(:,i,1));d1(i)=d;[~,d]=max(Mt(:,i,2));d2(i)=d;end
% d=d1-d2;
%figure, plot(1:numThetan,d,'r.-',1:numThetan,round(shift),'b.-');xlabel('\theta');ylabel('shift corresponding to error-free');legend('actual','analytical')
