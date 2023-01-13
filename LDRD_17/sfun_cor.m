function [f,g,XTM,shift]=sfun_COR(x,XTM,Ltol) 
global frame n_delta 
global N thetan numThetan dTau nTau 
global sinoS I0 
%%===== Reconstruction discrete objective
%%===== Ltol: intersection length matrix
%%===== f: sum_i ||e^T(Ltol_i.*I)e-M_i||^2, i=1..theta
%%===== x(1:n_delta): off center for the initial reference projection;

theta=thetan*pi/180;
if(n_delta==2*numThetan)
    %%=============================== COR for each projection
    shift=(cos(theta)-1).*x(1:2:n_delta)'+sin(theta).*x(2:2:n_delta)';
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
for i = 1:numThetan
    delay=shift(i);
    %%%=============================================
    % interpolate_rate=2;
    % H=fft(interp1(1:nTau+1,XTM(i,:),1:1/interpolate_rate:nTau+1,'nearest'));
    % u=[0:1/interpolate_rate:floor((nTau+1)/2)-1 floor(-(nTau+1)/2):1/interpolate_rate:-0.5] / (2*nTau+1);
    % ubar=exp(-j*2*pi.*(u*delay*interpolate_rate*2));
    % H1=H.*ubar;
    % align_temp=real(ifft(H1));
    % alignedSignal(i,:)=align_temp(1:interpolate_rate:end);
    % Dalign_temp1=real(ifft(H.*ubar.*(-2*pi*j*u*4)));
    % DalignedSignal(i,:)=Dalign_temp1(1:interpolate_rate:end);
    %%%=============================================
    H=fft(XTM(i,:));
    u=[0:floor((nTau+1)/2)-1 floor(-(nTau+1)/2):-1] / (nTau+1);
    ubar=exp(-j*2*pi.*(u*delay));
    alignedSignal(i,:)=abs(real(ifft(H.*ubar)));
    DalignedSignal(i,:)=sign(real(ifft(H.*ubar))).*real(ifft(H.*ubar.*(-2*pi*j*u)));
    %%%=============================================
    % H=fftshift(fft(XTM(i,:)));
    % u=[floor(-(nTau+1)/2):floor((nTau+1)/2)-1]/(nTau+1);
    % ubar=exp(-j*2*pi.*(u*delay));
    % alignedSignal1(i,:)=real(ifft(fftshift(H.*ubar)));
    % DalignedSignal(i,:)=real(ifft(fftshift(H.*ubar.*(-2*pi*j*u))));
    % % %%%=============================================
    % subplot(1,2,1),plot(1:nTau+1,alignedSignal(i,:),'r.-');hold on; plot(1:nTau+1,XTM(i,:),'b.-');hold off;title(num2str(delay));
    % subplot(1,2,2),plot(1:nTau+1,alignedSignal(i,:),'r.-');hold on; plot(1:nTau+1,alignedSignal1(i,:),'b.-');hold off;title(num2str(delay));pause;%,1:nTau+1,alignedSignal1(i,:),'g.-'
end
%%------------------------------------------------------
Daligned=zeros(numThetan*(nTau+1),n_delta);
for i=1:n_delta/2
    Daligned(i:numThetan:end,2*i-1)=DalignedSignal(i,:);
    Daligned(i:numThetan:end,2*i)=DalignedSignal(i,:);
end

Daligned=Daligned.*repmat(Ddelta,[1,nTau+1])';
% save('sfun_test.mat','Daligned','Ddelta','alignedSignal','shift');
% figure, subplot(1,2,1),imagesc(alignedSignal');subplot(1,2,2);imagesc(reshape(x(n_delta+1:end),N,N));
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
penalty=0;
if(penalty & n_delta==numThetan*2)
    % Tik=delsq(numgrid('S',N(1)+2)); 
    [~,~,Tik]=laplacian(n_delta/2,{'DD'});
    L1_norm=0;
    L2_norm=1;
    if(L2_norm)
        lambda=1e2;
        Reg=Tik*shift';
        f=f+lambda*(sum(Reg.^2));
        g_temp=[(cos(theta')-1).*(Tik'*Tik*shift'),sin(theta').*(Tik'*Tik*shift')]';
        g=g+lambda*2*[g_temp(:);zeros(N^2,1)];
        % Reg=Tik*x(n_delta+1:end);
        % f=f+lambda*(norm(Reg))^2;
        % g=g+lambda*2*[zeros(n_delta,1);Tik'*Reg];
    elseif(L1_norm)
        Reg=sum(abs(W(:)));
        f=f+lambda*Reg;
        g=g+lambda;
    end
end


