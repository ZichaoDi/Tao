%% Analytical formula to find the shift of each projection due to center of rotation transformation;
% c=15;
% delta(1)= c*Tau/1748;
% delta(2)= c*Tau/1748;

delta_d=0; % off center for the initial reference projection;
theta=thetan*pi/180;
step=-3:1e-6:3;
err_step=zeros(length(step),1);
shift=zeros(length(step),numThetan);
eps1=zeros(numThetan,1);
for k=1:length(step)
    k_step=step(k); 
    delta=delta0+[k_step,k_step].*dz;
    Num=(DetKnot0(1,1)-SourceKnot0(1,1))*(cos(theta)*delta(2)-sin(theta)*delta(1)+2*cos(theta).*sin(theta)*delta(1)+(sin(theta).^2-cos(theta).^2)*delta(2))...
    + (DetKnot0(1,2)-SourceKnot0(1,2))*(-sin(theta)*delta(2)-cos(theta)*delta(1)+2*cos(theta).*sin(theta)*delta(2)+(cos(theta).^2-sin(theta).^2)*delta(1));
    Det=norm(DetKnot0(1,:)-SourceKnot0(1,:));
    shift(k,:)=round(Num./Det/dTau)+delta_d;

    Mt=-log(DisR./I0);
    alignedSignal=zeros(nTau+1,numThetan);
    tt=Mt(:,:,2);
    for i = 1: numThetan
        delay=shift(i);
        if(delay>=0)
            alignedSignal(delay+1:end,i)=tt(1:end-delay,i);
            alignedSignal(1:delay,i)=tt(end-delay+1:end,i);
        else
            alignedSignal(1:end+delay,i)=tt(-delay+1:end,i);
            alignedSignal(end+delay+1:end,i)=tt(1:-delay,i);
        end
    end
    % err=norm(alignedSignal(:,3)-Mt(:,3,1));
    err_step(k)=err;
end
for i = 1:numThetan
    v_temp=unique(shift(:,i));
    eps1(i)=length(find(shift(:,i)==v_temp(round(length(v_temp)/2))));
end
% fprintf('Correcting Error = %d\n',err);
%  figure, 
%  subplot(3,2,1);imagesc(tt);subplot(3,2,2);imagesc(iradon(tt,thetan,'linear','shepp-logan',N));
%  subplot(3,2,3);imagesc(alignedSignal);subplot(3,2,4);imagesc(iradon(alignedSignal,thetan,'linear','shepp-logan',N));
%  subplot(3,2,5);imagesc(Mt(:,:,1));subplot(3,2,6);imagesc(iradon(Mt(:,:,1),thetan));


% d1=[];d2=[];for i=1:numThetan;[~,d]=max(Mt(:,i,1));d1(i)=d;[~,d]=max(Mt(:,i,2));d2(i)=d;end
% d=d1-d2;
%figure, plot(1:numThetan,d,'r.-',1:numThetan,round(shift),'b.-');xlabel('\theta');ylabel('shift corresponding to error-free');legend('actual','analytical')
