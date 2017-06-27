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
    shift(k,:)=Num./Det/dTau+1/2*delta_d0/dTau;
    Mt=-log(DisR./I0);%DisR;%
    alignedDiscrete=zeros(nTau+1,numThetan);
    LR2=zeros(nTau+1,numThetan);
    alignedContinuous=alignedDiscrete;
    tt=Mt(:,:,2);
    for i = 1:numThetan
        % delay=floor(mod(shift(k,i),nTau+1));
        % delay=double(mod(int32(shift(k,i)),nTau+1));
        delay=shift(k,i);
        %%---------------- Continuous: Gaussian Convolution
        sigma=0.005;%1.5/2.355;
        interpolate_rate=0.1;
        %%=============================================
        x_coarse=[-nTau-1:nTau+1];
        x_fine=[-nTau-1:interpolate_rate:nTau+1];
        tt_coarse=[zeros(nTau+2,1);tt(:,i)];
        tt_fine=interp1(x_coarse,tt_coarse,x_fine)';
        %%=============================================
        G=exp(-(x_fine'-delay).^2./(2*sigma^2));
        G=G/max(G);
        G_temp=fft(G).*fft(tt_fine);
        aligned_temp=ifft(fft(G).*fft(tt_fine));
        aligned_temp=aligned_temp(1:1/interpolate_rate:end);
        alignedContinuous(:,i)=aligned_temp(1:nTau+1);
        %%=========================================================
        uu=interp1([1:nTau+1],tt(:,i),[1:0.5:nTau+1]);
        H=fft(uu);
        u=[0:0.5:floor((nTau+1)/2)-1 floor(-(nTau+1)/2):0.5:-0.5 ]' / (2*nTau+1);
        u=exp(-sqrt(-1)*2*pi.*(uu*delay*2));
        % u((nTau+1)/2+1)=real(u((nTau+1)/2+1));
        H1=H.*u;
        temp=real(ifft(H1));
        LR2(:,i)=temp(1:2:end);
        % H=fft(tt(:,i));
        % u=[0:floor((nTau+1)/2)-1 floor(-(nTau+1)/2):-1]' / (nTau+1);
        % u=exp(-sqrt(-1)*2*pi.*(u*delay));
        % % u((nTau+1)/2+1)=real(u((nTau+1)/2+1));
        % H1=H.*u;
        % LR2(:,i)=real(ifft(H1));
        % figure,plot(LR2(:,i),'r.-');hold on; plot(tt(:,i),'b.-');title(num2str(delay));pause(1);
        %%=========================================================
        % x_coarse=[1:nTau+1];
        % x_fine=[1:interpolate_rate:nTau+1];
        % tt_coarse=[tt(:,i)];
        % tt_fine=interp1(x_coarse,tt_coarse,x_fine)';
        % G=exp(-(x_fine'-delay).^2./(2*sigma^2));
        % aligned_temp=ifft(fft(G).*fft(tt_fine));
        % figure,plot(aligned_temp,'r.-');pause(1);
        % alignedContinuous(:,i)=aligned_temp(1:1/interpolate_rate:end);
        %%=========================================================
        % G=exp(-([-nTau-1:nTau+1]'-delay).^2./(2*sigma^2));
        % G=G./max(G);
        % aligned_temp=ifft(fft(G).*fft([0*ones(nTau+2,1);tt(:,i)]));
        % aligned_temp=aligned_temp./max(aligned_temp)*max(tt(:,i));
        % if(abs(delay)>(nTau+1)/2)
        %     alignedContinuous(:,i)=aligned_temp(nTau+2:end-1);
        % else
        %     alignedContinuous(:,i)=aligned_temp(1:nTau+1);
        % end
        % if(delay>=0)
        %     alignedContinuous(:,i)=aligned_temp(nTau+2:end-1);
        % else
        %    alignedContinuous(:,i)=aligned_temp(1:nTau+1);
        % end
        %%---------------- Discrete: Permutation Matrix
        if(delay>=0)
            alignedDiscrete(delay+1:end,i)=tt(1:end-delay,i);
            alignedDiscrete(1:delay,i)=tt(end-delay+1:end,i);
        else
            alignedDiscrete(1:end+delay,i)=tt(-delay+1:end,i);
            alignedDiscrete(end+delay+1:end,i)=tt(1:-delay,i);
        end
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
subplot(nrow,2,[1,2]);imagesc(tt); title('Observed Sinogram')
% axis xy image;
% subplot(nrow,2,2);imagesc(iradon(fliplr(tt),thetan,'linear','shepp-logan',N)'); 
axis xy image;
subplot(nrow,2,[3,4]);imagesc(alignedContinuous);title('Continuous Corrected')
% subplot(nrow,2,4);imagesc(iradon(fliplr(alignedContinuous),thetan)');%,'linear','shepp-logan',N)'); 
axis xy image;
subplot(nrow,2,[5,6]);imagesc(alignedDiscrete);title('Discrete Corrected')
% subplot(nrow,2,6);imagesc(iradon(fliplr(alignedDiscrete),thetan)');%,'linear','shepp-logan',N)'); 
axis xy image;
subplot(nrow,2,7);imagesc(Mt(:,:,1));title('Sinogram with COR(0,0)')
% axis xy image;
subplot(nrow,2,8);imagesc(W);axis xy image;hold on;
plot(delta0(:,1)/dz(1)+N/2,delta0(:,2)/dz(2)+N/2,'r.-'); 
% for i=1:numThetan
%     text(delta0(i,1)/dz(1)+N/2,delta0(i,2)/dz(2)+N/2,num2str(i));
% end
% axis xy image;


% d1=[];d2=[];for i=1:numThetan;[~,d]=max(Mt(:,i,1));d1(i)=d;[~,d]=max(Mt(:,i,2));d2(i)=d;end
% d=d1-d2;
%figure, plot(1:numThetan,d,'r.-',1:numThetan,round(shift),'b.-');xlabel('\theta');ylabel('shift corresponding to error-free');legend('actual','analytical')
