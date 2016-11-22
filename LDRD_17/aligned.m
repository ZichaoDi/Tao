function alignedSignal=aligned(delta,Mt,Det0,Source0,thetan,nTau,numThetan,dTau)

delta_d=0; % off center for the initial reference projection;
theta=thetan*pi/180;
Num=(Det0(1)-Source0(1))*(cos(theta)*delta(2)-sin(theta)*delta(1)+2*cos(theta).*sin(theta)*delta(1)+(sin(theta).^2-cos(theta).^2)*delta(2))...
+ (Det0(2)-Source0(2))*(-sin(theta)*delta(2)-cos(theta)*delta(1)+2*cos(theta).*sin(theta)*delta(2)+(cos(theta).^2-sin(theta).^2)*delta(1));
Det=norm(Det0-Source0);
shift=mod(round(Num./Det/dTau)+delta_d,nTau+1);

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

