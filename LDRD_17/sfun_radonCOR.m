function [f,g,r]=sfun_radonCOR(x,XTM,Ltol) 
global frame  
global thetan numThetan nTau dTau DetKnot0 SourceKnot0 
%%===== Reconstruction discrete objective
%%===== Ltol: intersection length matrix
%%===== f: sum_i ||e^T(Ltol_i.*I)e-M_i||^2, i=1..theta
delta_d=0; % off center for the initial reference projection;
delta=x(1:2);
theta=thetan*pi/180;

Det=norm(DetKnot0(1,:)-SourceKnot0(1,:));
Num=(DetKnot0(1,1)-SourceKnot0(1,1))*(cos(theta).*delta(2)'-sin(theta).*delta(1)'+2*cos(theta).*sin(theta).*delta(1)'+(sin(theta).^2-cos(theta).^2).*delta(2)')...
+ (DetKnot0(1,2)-SourceKnot0(1,2))*(-sin(theta).*delta(2)'-cos(theta).*delta(1)'+2*cos(theta).*sin(theta).*delta(2)'+(cos(theta).^2-sin(theta).^2).*delta(1)');
shift=Num./Det/dTau+delta_d;

Ddelta=[(DetKnot0(1,1)-SourceKnot0(1,1))*(-sin(theta)+sin(2*theta))+(DetKnot0(1,2)-SourceKnot0(1,2))*(cos(2*theta)-cos(theta));...
(DetKnot0(1,1)-SourceKnot0(1,1))*(cos(theta)-cos(2*theta))+(DetKnot0(1,2)-SourceKnot0(1,2))*(sin(2*theta)-sin(theta))];
Ddelta=Ddelta./(Det*dTau);

alignedSignal=zeros(numThetan,nTau+1);
DalignedSignal=zeros(numThetan,nTau+1);
for i = 1:numThetan
    delay=shift(i);
    sigma=1;%length(find(tt(:,i)>0))/2.355;
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

Daligned=repmat(DalignedSignal(:),[1,2]).*repmat(Ddelta,[1,nTau+1])';
XTM=alignedSignal;

if(strcmp(frame,'EM'))
    thres=1;
    Mt=XTM(:)+thres;
    Rdis=Ltol*x(3:end)+thres;
    f=sum(-log(Rdis).*Mt+Rdis);
    g= Ltol'*(-Mt./Rdis+1);
    g_pert=-log(Rdis)'*(Daligned);
elseif(strcmp(frame,'LS'))
    r=Ltol*x(3:end)-XTM(:);
    g=Ltol'*r;
    f=1/2*sum(r.^2,1);
    g_pert=-r'*Daligned;
end
g=[g_pert';g];




