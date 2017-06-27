function [f,g,XTM]=sfun_cor_dr(x,xstar,XTM,Ltol) 
global frame N_delta 
global thetan numThetan dTau nTau 
global sinoS N
%%===== Reconstruction discrete objective
%%===== Ltol: intersection length matrix
%%===== f: sum_i ||e^T(Ltol_i.*I)e-M_i||^2, i=1..theta
%%===== x(1:N_delta): off center for the initial reference projection;

shift=x;
alignedSignal=zeros(numThetan,nTau+1);
DalignedSignal=zeros(numThetan,nTau+1);
sigma=0.05;%1.5/2.355;
for i = 1:numThetan
    % delay=shift(i);
    delay=double(mod(int32(shift(i)),nTau+1));

    G=exp(-([-nTau-1:nTau+1]'-delay).^2./(2*sigma^2));
    dG=(([-nTau-1:nTau+1]'-delay)./(sigma^2)).*exp(-([-nTau-1:nTau+1]'-delay).^2./(2*sigma^2));
    aligned_temp=ifft(fft(G).*fft([zeros(nTau+2,1);XTM(i,:)']));
    Daligned=ifft(fft(dG).*fft([zeros(nTau+2,1);XTM(i,:)']));
    if(abs(delay)>(nTau+1)/2)
        alignedSignal(i,:)=aligned_temp(nTau+2:end-1);
        DalignedSignal(i,:)=Daligned(nTau+2:end-1);
    else
        alignedSignal(i,:)=aligned_temp(1:nTau+1);
        DalignedSignal(i,:)=Daligned(1:nTau+1);
    end
end

%%------------------------------------------------------
Daligned=zeros(numThetan*(nTau+1),N_delta);
for i=1:N_delta
    Daligned(i:numThetan:end,i)=DalignedSignal(i,:);
end
% Daligned1=repmat(DalignedSignal(:),[1,N_delta]);
% save('sfun_test.mat','Daligned','alignedSignal','shift');
% figure, subplot(1,2,1),imagesc(alignedSignal');subplot(1,2,2);imagesc(reshape(x(N_delta+1:end),N,N));
XTM=alignedSignal;

if(strcmp(frame,'EM'))
    thres=1;
    Mt=XTM(:)+thres;
    Rdis=Ltol*xstar+thres;
    f=sum(-log(Rdis).*Mt+Rdis);
    g=-Daligned'*log(Rdis);
elseif(strcmp(frame,'LS'))
    r=Ltol*x(N_delta+1:end)-XTM(:);
    f=1/2*sum(r.^2,1);
    g=-Daligned'*r;
end
