function [f,g,XTM,f_sub]=sfun_shift_mc(x,XTM,Ltol) 
global frame N_delta  
global thetan numThetan dTau nTau 
global lambda sinoS N NumElement alpha
global reg_str
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
Tik=delsq(numgrid('S',N+2)); 
for ele=1:NumElement
    if(strcmp(frame,'EM'))
        thres=1;
        Mt=XTM(:,ele)+thres;
        Rdis=Ltol*xc(:,ele)+thres;
        if(strcmp(reg_str,'TV'))
            epsilon=1e-4;
            [~,~,Tik]=laplacian(N,{'DD'});
            x = reshape(xc(:,ele),N,N);
            f2_v=sqrt((Tik*x).^2 + (Tik*x').^2+ epsilon);
            f2=sum(sum(f2_v));
            f_sub(ele)=sum(-log(Rdis).*Mt+Rdis);
            f=f+alpha(ele)*sum(-log(Rdis).*Mt+Rdis)+lambda*f2;
            g_reg= Tik'*Tik*x./f2_v;
            g= [g; alpha(ele)*Ltol'*(-Mt./Rdis+1)+lambda*g_reg(:)];
        elseif(strcmp(reg_str,'L2'))
            %%============L2
            reg=Tik*xc(:,ele);
            f=f+alpha(ele)*sum(-log(Rdis).*Mt+Rdis)+lambda*sum(reg.^2);
            g= [g; alpha(ele)*Ltol'*(-Mt./Rdis+1)+lambda*2*Tik'*Tik*xc(:,ele)];
        end
        g_pert=g_pert-alpha(ele)*log(Rdis)'*(Daligned(:,:,ele))';
    elseif(strcmp(frame,'LS'))
        r=Ltol*xc(:,ele)-XTM(:,ele);
        if(strcmp(reg_str,'TV'))
            %%============ fake TV
            epsilon=1e-4;
            [~,~,Tik]=laplacian(N,{'DD'});
            x = reshape(xc(:,ele),N,N);
            f2_v=sqrt((Tik*x).^2 + (Tik*x').^2+ epsilon);
            f2=sum(sum(f2_v));
            f=f+1/2*alpha(ele)*r'*r+lambda*f2;
            g_reg= Tik'*Tik*x./f2_v;
            g=[g;alpha(ele)*Ltol'*r+lambda*g_reg(:)];
            g_pert=g_pert-alpha(ele)*r'*(Daligned(:,:,ele))';
        elseif(strcmp(reg_str,'L2'))
            %%============L2
            reg=Tik*xc(:,ele);
            f=f+1/2*alpha(ele)*sum(r.^2,1)+lambda*sum(reg.^2);
            g=[g;alpha(ele)*Ltol'*r+lambda*2*Tik'*Tik*xc(:,ele)];
            g_pert=g_pert-alpha(ele)*r'*(Daligned(:,:,ele))';
        end
    end
end
g=[g_pert';g];
