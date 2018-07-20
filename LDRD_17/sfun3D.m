function [f,g]=sfun3D(x,project,Ltol) 
global totN frame numThetan nTau nslice N NumElement
global alpha
xc=reshape(x,N^2,NumElement,nslice);
f=0;
%%%=======================================================
g=zeros(N^2,NumElement,nslice);
XTM=reshape(project,numThetan*(nTau+1),nslice,NumElement);
for ele=1:NumElement
    if(strcmp(frame,'EM'))
        thres=1;
        Mt=XTM(:,:,ele)+thres;
        Rdis=Ltol*squeeze(xc(:,ele,:))+thres;
        f=f+alpha(:,ele)'*sum(-log(Rdis).*Mt+Rdis,1)';
        g(:,ele,:)= Ltol'*(-Mt./Rdis+1).*repmat(alpha(:,ele)',N^2,1);
    elseif(strcmp(frame,'LS'))
        r=Ltol*squeeze(xc(:,ele,:))-XTM(:,:,ele);
        g(:,ele,:)=Ltol'*r.*repmat(alpha(:,ele)',N^2,1);
        f=f+1/2*sum(alpha(:,ele)'.*sum(r.^2,1),2);
    end
end
scale=1;
f=f./scale;
g=g(:)./scale;
