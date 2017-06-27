function [f,g,XTM]=sfun_COR_dr(x,XTM,Ltol) 
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
sigma=0.05;%1.5/2.355;
for i = 1:numThetan
    delay=shift(i);
    % % delay=double(mod(int32(shift(i)),nTau+1));
    % x_coarse=[-nTau-1:nTau+1];
    % interpolate_rate=1;
    % x_fine=[-nTau-1:interpolate_rate:nTau+1];
    % tt_coarse=[zeros(nTau+2,1);XTM(i,:)'];
    % tt_fine=interp1(x_coarse,tt_coarse,x_fine)';
    % G=exp(-(x_fine'-delay).^2./(2*sigma^2));
    % dG=((x_fine'-delay)./(sigma^2)).*exp(-(x_fine'-delay).^2./(2*sigma^2));
    % % figure, subplot(1,2,1),plot(G,'r.-');subplot(1,2,2);plot(dG,'b.-'); pause;
    % % G=1/(sigma*sqrt(2*pi))*exp(-([-nTau-1:nTau+1]'-delay).^2./(2*sigma^2));
    % % dG=1/(sigma*sqrt(2*pi))*(([-nTau-1:nTau+1]'-delay)./(sigma^2)).*exp(-([-nTau-1:nTau+1]'-delay).^2./(2*sigma^2));
    % aligned_temp=ifft(fft(G).*fft(tt_fine));
    % aligned_temp=aligned_temp(1:1/interpolate_rate:end);
    % Daligned=ifft(fft(dG).*fft(tt_fine));
    % % figure,plot(Daligned,'r.-')
    % %%%=============================================
    % Daligned=Daligned(1:1/interpolate_rate:end);
    % alignedSignal(i,:)=aligned_temp(1:nTau+1);
    % DalignedSignal(i,:)=Daligned(1:nTau+1);
    % %%%=============================================
        H=fft(XTM(i,:));
        u=[0:floor((nTau+1)/2)-1 floor(-(nTau+1)/2):-1] / (nTau+1);
        ubar=exp(-j*2*pi.*(u*delay));
        % u((nTau+1)/2+1)=real(u((nTau+1)/2+1));
        H1=H.*ubar;
        alignedSignal(i,:)=real(ifft(H1));
        DalignedSignal(i,:)=real(ifft(H.*ubar.*(-2*pi*j*u)));
end

%%------------------------------------------------------
Daligned=zeros(numThetan*(nTau+1),N_delta);
for i=1:N_delta
    Daligned(i:numThetan:end,i)=DalignedSignal(i,:);
end
% Daligned1=repmat(DalignedSignal(:),[1,N_delta]);
% save('sfun_test.mat','Daligned','alignedSignal','shift');
% figure, subplot(1,2,1),imagesc(alignedSignal');subplot(1,2,2);imagesc(reshape(x(N_delta+1:end),N,N));
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
    g=Ltol'*r;
    f=1/2*sum(r.^2,1);
    g_pert=-r'*Daligned;
end
g=[g_pert';g];
penalty=0;
if(penalty & N_delta==numThetan)
    % Tik=delsq(numgrid('S',N(1)+2)); 
    [~,~,Tik]=laplacian(N_delta,{'DD'});
    L1_norm=0;
    L2_norm=1;
    if(L2_norm)
        lambda=1e1;
        Reg=Tik*x(1:N_delta);
        f=f+lambda*(norm(Reg))^2;
        g=g+lambda*2*[Tik'*Reg;zeros(N^2,1)];
        % Reg=Tik*x(n_delta+1:end);
        % f=f+lambda*(norm(Reg))^2;
        % g=g+lambda*2*[zeros(n_delta,1);Tik'*Reg];
    elseif(L1_norm)
        Reg=sum(abs(W(:)));
        f=f+lambda*Reg;
        g=g+lambda;
    end
end


