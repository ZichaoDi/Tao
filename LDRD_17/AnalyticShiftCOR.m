%% Analytical formula to find the shift of each projection due to center of rotation transformation;
% c=15;
% delta(1)= c*Tau/1748;
% delta(2)= c*Tau/1748;
global Ddelta_COR

delta_d=0; % off center for the initial reference projection;
theta=thetan*pi/180;
step=[0];%-10:1e-1:10;
eps1=zeros(numThetan,1);

Det=norm(DetKnot0(1,:)-SourceKnot0(1,:));

eps2_Det=Det*dTau;
eps2_Num=...
[2*(DetKnot0(1,1)-SourceKnot0(1,1))*(-sin(theta)+sin(2*theta))+2*(DetKnot0(1,2)-SourceKnot0(1,2))*(cos(2*theta)-cos(theta)); ...
2*(DetKnot0(1,1)-SourceKnot0(1,1))*(cos(theta)-cos(2*theta))+2*(DetKnot0(1,2)-SourceKnot0(1,2))*(sin(2*theta)-sin(theta))];
eps2=eps2_Det./eps2_Num./dz(1);
eps2(isinf(eps2))=0;

Ddelta=[(DetKnot0(1,1)-SourceKnot0(1,1))*(-sin(theta)+sin(2*theta))+(DetKnot0(1,2)-SourceKnot0(1,2))*(cos(2*theta)-cos(theta));...
(DetKnot0(1,1)-SourceKnot0(1,1))*(cos(theta)-cos(2*theta))+(DetKnot0(1,2)-SourceKnot0(1,2))*(sin(2*theta)-sin(theta))];
Ddelta=Ddelta/eps2_Det;
Ddelta_COR=repmat(Ddelta,[1 nTau+1]);
k_outer=[1];
err_step=zeros(length(k_outer),length(step));
shift=zeros(length(step),numThetan);
for ko_ind=1:length(k_outer);

delta0=k_outer(ko_ind)*cr(2,:);
for k=1:length(step)
    k_step=step(k); 
    delta=repmat(delta0',[1,numThetan])+pert';%[repmat(cr(2,:),[floor(numThetan/2),1]);repmat(cr(3,:),[numThetan-floor(numThetan/2),1])]';
    % +k_step*Ddelta'.*abs(ceil(eps2))'.*repmat(dz,[numThetan,1]);
    Num=(DetKnot0(1,1)-SourceKnot0(1,1))*((cos(theta)-1).*delta(2,:)+delta(1,:).*sin(theta))+ (DetKnot0(1,2)-SourceKnot0(1,2))*((1-cos(theta)).*delta(1,:)+sin(theta).*delta(2,:));
    shift(k,:)=Num./Det/dTau+1/2*delta_d0/dTau;
    Mt=-log(DisR./I0);%DisR;%
    alignedDiscrete=zeros(nTau+1,numThetan);
    alignedContinuous=alignedDiscrete;
    tt=Mt(:,:,2);
    sigma=1.5/2.355;
    for i = 1:numThetan
        delay=floor(mod(shift(k,i),nTau+1));
        %%---------------- Continuous: Gaussian Convolution
        G=1/(sigma*sqrt(2*pi))*exp(-([-nTau-1:nTau+1]'-delay).^2./(2*sigma^2));
        aligned_temp=ifft(fft(G).*fft([0*ones(nTau+2,1);tt(:,i)]));
        if(abs(delay)>(nTau+1)/2)
            alignedContinuous(:,i)=aligned_temp(nTau+3:end);
        else
            alignedContinuous(:,i)=aligned_temp(1:nTau+1);
        end
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
% err(sig_i,:)=[norm(alignedContinuous-Mt(:,:,1)),norm(alignedDiscrete-Mt(:,:,1))];
% for i = 1:numThetan
%     v_temp=unique(shift(:,i));
%     eps1(i)=(step(2)-step(1))*length(find(shift(:,i)==v_temp(round(length(v_temp)/2))));
% end
% fprintf('Correcting Error = %d\n',err);
figure, 
subplot(3,2,1);imagesc(tt); title('Sinogram with True COR')
% axis xy image;
subplot(3,2,2);imagesc(iradon(fliplr(tt),thetan,'linear','shepp-logan',N)'); 
axis xy image;
subplot(3,2,3);imagesc(alignedContinuous);title('Corrected Sinogram')
% axis xy image;
subplot(3,2,4);imagesc(iradon(fliplr(alignedContinuous),thetan)');%,'linear','shepp-logan',N)'); 
axis xy image;
subplot(3,2,5);imagesc(Mt(:,:,1));title('Sinogram with COR(0,0)')
% axis xy image;
subplot(3,2,6);imagesc(W);axis xy image;hold on;
plot(delta(1,:)/dz(1)+N/2,delta(2,:)/dz(2)+N/2,'r.-'); 
% axis xy image;


% d1=[];d2=[];for i=1:numThetan;[~,d]=max(Mt(:,i,1));d1(i)=d;[~,d]=max(Mt(:,i,2));d2(i)=d;end
% d=d1-d2;
%figure, plot(1:numThetan,d,'r.-',1:numThetan,round(shift),'b.-');xlabel('\theta');ylabel('shift corresponding to error-free');legend('actual','analytical')
