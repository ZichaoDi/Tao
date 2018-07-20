function [f,g,XTM]=sfun_cor_dr(x,XTM,Ltol) 
global frame N_delta 
global thetan numThetan dTau nTau 
global sinoS N
%%===== Reconstruction discrete objective
%%===== Ltol: intersection length matrix
%%===== f: sum_i ||e^T(Ltol_i.*I)e-M_i||^2, i=1..theta
%%===== x(1:N_delta): off center for the initial reference projection;

shift=x(1:N_delta);
sigma=1.5/2.355;
scale=1/sqrt(2*pi)/sigma;
extNtau=nTau+1+20;
realInd=11:nTau+1+10;
padInd=setdiff([1:extNtau],realInd);
alignedSignal=zeros(numThetan,nTau+1);
DalignedSignal=zeros(numThetan,extNtau);
for i = 1:numThetan
    temp=zeros(extNtau,1);temp(realInd)=XTM(i,:);
    temp1=gaussfilt(1:extNtau,temp,5/2.355);
    temp(padInd)=temp1(padInd); 
    delay=shift(i);
    range=[0:ceil(extNtau/2)-1 ceil(-extNtau/2):-1]';%
    G=exp(-(range-delay).^2./(2*sigma^2));
    dG=((range-delay)./(sigma^2)).*exp(-(range-delay).^2./(2*sigma^2));
    aligned_temp=scale*real(ifft(fft(G).*(fft(temp))));
    alignedSignal(i,:)=aligned_temp(realInd);
    DalignedSignal(i,:)=scale*real(ifft(fft(dG).*(fft(temp))));
    % plot(1:extNtau,temp,'m-',realInd,alignedSignal(i,:),'r.-',realInd,XTM(i,:),'g.-',1:extNtau,aligned_temp,'b-',[11,11],[0,1],'k--',[nTau+11,nTau+11],[0,1],'k--')
    % legend('original padded','shifted','original','shifted padded')
    % title(num2str([i,shift(i)]));
    % pause;
end
%%------------------------------------------------------
Daligned=zeros(numThetan*(nTau+1),N_delta);
for i=1:N_delta
    Daligned(i:numThetan:end,i)=DalignedSignal(i,realInd);
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
    g_pert=-r'*Daligned;
end
g=[g_pert';g];
