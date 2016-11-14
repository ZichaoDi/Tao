function [f,g,r]=sfun_radonCOR(x,Mt,Ltol) 
global frame 
global thetan numThetan nTau dTau DetKnot0 SourceKnot0 
%%===== Reconstruction discrete objective
%%===== Ltol: intersection length matrix
%%===== f: sum_i ||e^T(Ltol_i.*I)e-M_i||^2, i=1..theta
MU=x(1:end);
% alignedSignal=alinged(x(1:2),reshape(Mt,nTau+1,numThetan),DetKnot0(1,:),SourceKnot0(1,:),thetan,nTau,numThetan,dTau);
alignedSignal=Mt;
if(strcmp(frame,'EM'))
    thres=1;
    Mt=alignedSignal(:)+thres;
    Rdis=Ltol*MU+thres;
    f=sum(-log(Rdis).*Mt+Rdis);
    g= Ltol'*(-Mt./Rdis+1);
elseif(strcmp(frame,'LS'))
    Mt=alignedSignal(:);
    r=Ltol*MU-Mt;
    g=Ltol'*r;
    f=1/2*sum(r.^2,1);
end


function alignedSignal=alinged(delta,Mt,Det0,Source0,thetan,nTau,numThetan,dTau)

delta_d=0; % off center for the initial reference projection;
theta=thetan*pi/180;
Num=(Det0(1)-Source0(1))*(cos(theta)*delta(2)-sin(theta)*delta(1)+2*cos(theta).*sin(theta)*delta(1)+(sin(theta).^2-cos(theta).^2)*delta(2))...
+ (Det0(2)-Source0(2))*(-sin(theta)*delta(2)-cos(theta)*delta(1)+2*cos(theta).*sin(theta)*delta(2)+(cos(theta).^2-sin(theta).^2)*delta(1));
Det=norm(Det0-Source0);
shift=round(Num./Det/dTau)+delta_d;

alignedSignal=zeros(nTau+1,numThetan);
for i = 1: numThetan
    delay=shift(i);
    if(delay>=0)
        alignedSignal(delay+1:end,i)=Mt(1:end-delay,i);
        alignedSignal(1:delay,i)=Mt(end-delay+1:end,i);
    else
        alignedSignal(1:end+delay,i)=Mt(-delay+1:end,i);
        alignedSignal(end+delay+1:end,i)=Mt(1:-delay,i);
    end
end

