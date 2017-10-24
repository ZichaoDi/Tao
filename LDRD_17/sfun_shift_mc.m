function [f,g,XTM]=sfun_shift_mc(x,XTM,Ltol) 
global frame N_delta 
global thetan numThetan dTau nTau 
global sinoS N NumElement alpha
%%===== Reconstruction discrete objective
%%===== Ltol: intersection length matrix
%%===== f: sum_i ||e^T(Ltol_i.*I)e-M_i||^2, i=1..theta
%%===== x(1:N_delta): off center for the initial reference projection;
%%===== XTM: multichannel sigogram

shift=x(1:N_delta);
xc=reshape(x(N_delta+1:end),N^2,NumElement);
alignedSignal=zeros(numThetan,nTau+1,NumElement);
DalignedSignal=zeros(numThetan,nTau+1,NumElement);
sigma=1.5/2.355;
scale=1/sqrt(2*pi)/sigma;
for ele=1:NumElement
for i = 1:numThetan
    delay=shift(i);
    range=[0:floor((nTau+1)/2)-1 floor(-(nTau+1)/2):-1];%
    G=exp(-(range-delay).^2./(2*sigma^2));
    dG=((range-delay)./(sigma^2)).*exp(-(range-delay).^2./(2*sigma^2));
    alignedSignal(i,:,ele)=scale*real(ifft(fft(G).*(fft(XTM(i,:,ele)))));
    DalignedSignal(i,:,ele)=scale*real(ifft(fft(dG).*(fft(XTM(i,:,ele)))));
end
end
%%------------------------------------------------------
Daligned=zeros(N_delta,numThetan*(nTau+1),NumElement);
for i=1:N_delta
    Daligned(i,i:numThetan:end,:)=DalignedSignal(i,:,:);
end
XTM=reshape(alignedSignal,numThetan*(nTau+1),NumElement);

f=0;
g_pert=zeros(1,N_delta);
g=[];
for ele=1:NumElement
if(strcmp(frame,'EM'))
    thres=1;
    Mt=XTM(:,ele)+thres;
    Rdis=Ltol*xc(:,ele)+thres;
    f=f+alpha(ele)*sum(-log(Rdis).*Mt+Rdis);
    g= [g; alpha(ele)*Ltol'*(-Mt./Rdis+1)];
    g_pert=g_pert-alpha(ele)*log(Rdis)'*(Daligned(:,:,ele))';
elseif(strcmp(frame,'LS'))
    r=Ltol*x(N_delta+1:end)-XTM(:);
    f=1/2*sum(r.^2,1);
    g=Ltol'*r;
    g_pert=-r'*Daligned;
end
end
g=[g_pert';g];
