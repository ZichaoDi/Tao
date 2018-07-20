function [f,g,y,Daligned]=sfun_shift_xy(x,project,Ltol) 
global totN frame numThetan nTau nslice N NumElement
global alpha
xc=reshape(x(2*numThetan+1:end),N^2,NumElement,nslice);
[y,Daligned]=GaussianShift2D(project,x(1:2*numThetan));
f=0;
g_pert=zeros(2*numThetan,1);
%%%=======================================================
g=zeros(N^2,NumElement,nslice);
XTM=reshape(y,numThetan*(nTau+1),nslice,NumElement);
for ele=1:NumElement
    if(strcmp(frame,'EM'))
        thres=1;
        Mt=XTM(:,:,ele)+thres;
        Rdis=Ltol*squeeze(xc(:,ele,:))+thres;
        f=f+alpha(:,ele)'*sum(-log(Rdis).*Mt+Rdis,1)';
        g(:,ele,:)= Ltol'*(-Mt./Rdis+1);
        gp_temp=log(Rdis).*repmat(alpha(:,ele)',numThetan*(nTau+1),1);
        g_pert(1:2:end)=g_pert(1:2:end)-(reshape(Daligned(1:2:end,:,:,:,ele),numThetan,numThetan*(nTau+1)*nslice))*gp_temp(:);
        g_pert(2:2:end)=g_pert(2:2:end)-(reshape(Daligned(2:2:end,:,:,:,ele),numThetan,numThetan*(nTau+1)*nslice))*gp_temp(:);
    elseif(strcmp(frame,'LS'))
        r=Ltol*squeeze(xc(:,ele,:))-XTM(:,:,ele);
        g(:,ele,:)=Ltol'*r;
        f=f+1/2*sum(sum(r.^2,1),2);
        g_pert(1:2:end)=g_pert(1:2:end)-(reshape(Daligned(1:2:end,:,:,:,ele),numThetan,numThetan*(nTau+1)*nslice))*r(:);
        g_pert(2:2:end)=g_pert(2:2:end)-(reshape(Daligned(2:2:end,:,:,:,ele),numThetan,numThetan*(nTau+1)*nslice))*r(:);
    end
end
g=[g_pert;g(:)];
