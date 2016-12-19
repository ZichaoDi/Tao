function [f,g,r]=sfun_radonCOR_sim(x,XTM,Ltol) 
global frame n_delta 
global thetan numThetan dTau nTau 
%%===== Reconstruction discrete objective
%%===== Ltol: intersection length matrix
%%===== f: sum_i ||e^T(Ltol_i.*I)e-M_i||^2, i=1..theta
%%===== x(1:3): off center for the initial reference projection;

theta=thetan*pi/180;

if(n_delta==2)
    shift=(cos(theta)-cos(2*theta))*x(1)+(sin(2*theta)-sin(theta))*x(2);
    Ddelta=[cos(theta)-cos(2*theta);sin(2*theta)-sin(theta)];
else
    shift=(cos(theta)-cos(2*theta))*x(1)+(sin(2*theta)-sin(theta))*x(2)+x(3);
    Ddelta=[cos(theta)-cos(2*theta);sin(2*theta)-sin(theta);ones(1,numThetan)];
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
% G=1/(sigma*sqrt(2*pi))*exp(-(repmat([1:nTau+1]',1,numThetan)-repmat(mod(shift,nTau+1),nTau+1,1)).^2./(2*sigma^2));
% dG=1/(sigma*sqrt(2*pi))*((repmat([1:nTau+1]',1,numThetan)-repmat(mod(shift,nTau+1),nTau+1,1))./(sigma^2)).*exp(-(repmat([1:nTau+1]',1,numThetan)-repmat(mod(shift,nTau+1),nTau+1,1)).^2./(2*sigma^2));
% alignedSignal1=(ifft(fft(G).*fft(XTM')))';
% DalignedSignal1=(ifft(fft(dG).*fft(XTM')))';
% norm(alignedSignal-alignedSignal1)
% norm(XTM-alignedSignal1)

%%------------------------------------------------------
Daligned=repmat(DalignedSignal(:),[1,n_delta]).*repmat(Ddelta,[1,nTau+1])';
% save('sfun_test.mat','alignedSignal','shift');
% figure, imagesc(iradon(alignedSignal',thetan));
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


