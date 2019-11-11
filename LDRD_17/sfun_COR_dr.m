function [f,g,XTM]=sfun_COR_dr(x,XTM,Ltol) 
global shift
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
for i = 1:numThetan
    delay=shift(i);
    % %%===========================================
    % interpolate_rate=2;
    % H=fft(interp1(1:nTau+1,XTM(i,:),1:1/interpolate_rate:nTau+1,'nearest'));
    % u=[0:1/interpolate_rate:floor((nTau+1)/2)-1 floor(-(nTau+1)/2):1/interpolate_rate:-0.5] / (interpolate_rate*nTau+1);
    % ubar=exp(-j*2*pi.*(u*delay*interpolate_rate*2));
    % H1=H.*ubar;
    % align_temp=real(ifft(H1));
    % alignedSignal(i,:)=align_temp(1:interpolate_rate:end);
    % Dalign_temp1=real(ifft(H.*ubar.*(-2*pi*j*u*2*interpolate_rate)));
    % DalignedSignal(i,:)=Dalign_temp1(1:interpolate_rate:end);
    % %%%=============================================
    H=fft(XTM(i,:));
    u=[0:floor((nTau+1)/2)-1 floor(-(nTau+1)/2):-1] / (nTau+1);
    ubar=exp(-j*2*pi.*(u*(delay)));
    alignedSignal(i,:)=real(ifft(H.*ubar));
    DalignedSignal(i,:)=real(ifft(H.*ubar.*(-2*pi*j*u)));
end

%%------------------------------------------------------
Daligned=zeros(numThetan*(nTau+1),N_delta);
for i=1:N_delta
    Daligned(i:numThetan:end,i)=DalignedSignal(i,:);
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
        lambda=1e5;
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


