function [f,g,XTM]=sfun_cor_dr(x,XTM,Ltol) 
global frame N_delta 
global thetan numThetan dTau nTau 
global sinoS N
%%===== Reconstruction discrete objective
%%===== Ltol: intersection length matrix
%%===== f: sum_i ||e^T(Ltol_i.*I)e-M_i||^2, i=1..theta
%%===== x(1:N_delta): off center for the initial reference projection;

shift=x(1:N_delta);
alignedSignal=zeros(numThetan,nTau+1);
DalignedSignal=zeros(numThetan,nTau+1);
sigma=0.001/2.355;
for i = 1:numThetan
    delay=shift(i);
    range=[0:floor((nTau+1)/2)-1 floor(-(nTau+1)/2):-1]';%
    G=exp(-(range-delay).^2./(2*sigma^2));
    dG=((range-delay)./(sigma^2)).*exp(-(range-delay).^2./(2*sigma^2));
    alignedSignal(i,:)=1*real(ifft(fft(G).*(fft(XTM(i,:)'))));
    DalignedSignal(i,:)=1*real(ifft(fft(dG).*(fft(XTM(i,:)'))));
    % subplot(1,2,1)
    % plot(1:nTau+1,alignedSignal(i,:),'r.-',1:nTau+1,sinoS(:,i),'b.-')
    % title(num2str(delay));
    % subplot(1,2,2)
    % plot(G,'r.-');
    % pause;
end
%%------------------------------------------------------
Daligned=zeros(numThetan*(nTau+1),N_delta);
for i=1:N_delta
    Daligned(i:numThetan:end,i)=DalignedSignal(i,:);
end
XTM=alignedSignal;

if(strcmp(frame,'EM'))
    thres=1;
    Mt=XTM(:)+thres;
    Rdis=Ltol*x(N_delta+1:end)+thres;
    f=sum(-log(Rdis).*Mt+Rdis);
    g= Ltol'*(-Mt./Rdis+1);
    g_pert=-log(Rdis)'*(Daligned);
elseif(strcmp(frame,'LS'))
    r=Ltol*x(N_delta+1:end)-XTM(:);
    f=1/2*sum(r.^2,1);
    g=Ltol'*r;
    g_pert=-Daligned'*r;
end
g=[g_pert';g];
