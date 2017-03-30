function [f,g,XTM]=sfun_COR(x,XTM,Ltol) 
global frame n_delta 
global thetan numThetan dTau nTau 
global sinoS
%%===== Reconstruction discrete objective
%%===== Ltol: intersection length matrix
%%===== f: sum_i ||e^T(Ltol_i.*I)e-M_i||^2, i=1..theta
%%===== x(1:3): off center for the initial reference projection;

theta=thetan*pi/180;
if(n_delta==2*numThetan)
    %%=============================== COR for each projection
    delta=reshape(x(1:numThetan*2),2,numThetan);
    shift=(cos(theta)-1).*delta(1,:)+sin(theta).*delta(2,:);
    i1=[1:2*numThetan];
    i2=reshape(repmat((1:numThetan),[2,1]),2*numThetan,1);
    i3=reshape([cos(theta)-1;sin(theta)],2*numThetan,1);
    Ddelta=sparse(i1,i2,i3); 
    %%========================== same COR for half of the projections
elseif(n_delta==4)
    delta(1,:)=[repmat(x(1),[1,floor(numThetan/2)]) repmat(x(3),[1,numThetan-floor(numThetan/2)])];
    delta(2,:)=[repmat(x(2),[1,floor(numThetan/2)]) repmat(x(4),[1,numThetan-floor(numThetan/2)])];
    shift=(cos(theta)-1).*delta(1,:)+sin(theta).*delta(2,:);
    Ddelta=zeros(n_delta,numThetan);
    Ddelta(1,1:floor(numThetan/2))=cos(theta(1:floor(numThetan/2)))-1;
    Ddelta(2,1:floor(numThetan/2))=sin(theta(1:floor(numThetan/2)));
    Ddelta(3,floor(numThetan/2)+1:numThetan)=cos(theta(floor(numThetan/2)+1:numThetan))-1;
    Ddelta(4,floor(numThetan/2)+1:numThetan)=sin(theta(floor(numThetan/2)+1:numThetan));
elseif(n_delta==2)
    shift=(cos(theta)-1)*x(1)+sin(theta)*x(2);
    Ddelta=[cos(theta)-1;sin(theta)];
end

alignedSignal=zeros(numThetan,nTau+1);
DalignedSignal=zeros(numThetan,nTau+1);
sigma=1.5/2.355;
for i = 1:numThetan
    delay=shift(i);
    G=1/(sigma*sqrt(2*pi))*exp(-([-nTau-1:nTau+1]'-delay).^2./(2*sigma^2));
    dG=1/(sigma*sqrt(2*pi))*(([-nTau-1:nTau+1]'-delay)./(sigma^2)).*exp(-([-nTau-1:nTau+1]'-delay).^2./(2*sigma^2));
    aligned_temp=ifft(fft(G).*fft([zeros(nTau+2,1);XTM(i,:)']));
    Daligned=ifft(fft(dG).*fft([zeros(nTau+2,1);XTM(i,:)']));
    if(abs(delay)>(nTau+1)/2)
        alignedSignal(i,:)=aligned_temp(nTau+3:end);
        DalignedSignal(i,:)=Daligned(nTau+3:end);
    else
        alignedSignal(i,:)=aligned_temp(1:nTau+1);
        DalignedSignal(i,:)=Daligned(1:nTau+1);
    end
end

%%------------------------------------------------------
Daligned=repmat(DalignedSignal(:),[1,n_delta]).*repmat(Ddelta,[1,nTau+1])';
% save('sfun_test.mat','Daligned','Ddelta','alignedSignal','shift');
% figure, subplot(1,2,1),imagesc(alignedSignal');subplot(1,2,2);imagesc(iradon(alignedSignal',thetan));
XTM=alignedSignal;

if(strcmp(frame,'EM'))
    thres=1;
    Mt=XTM(:)+thres;
    Rdis=Ltol*x(n_delta+1:end)+thres;
    f=sum(-log(Rdis).*Mt+Rdis);
    g= Ltol'*(-Mt./Rdis+1);
    g_pert=-log(Rdis)'*(Daligned);
elseif(strcmp(frame,'LS'))
    r=Ltol*x(n_delta+1:end)-XTM(:);
    g=Ltol'*r;
    f=1/2*sum(r.^2,1);
    g_pert=-r'*Daligned;
end
g=[g_pert';g];


