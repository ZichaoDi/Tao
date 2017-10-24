function [f,g]=sfun_shift_xy(x,project,Ltol) 
global totN frame numThetan nTau nslice N NumElement
global alpha
shift=reshape(x(1:2*numThetan),2,numThetan);
xc=reshape(x(2*numThetan+1:end),N^2,NumElement,nslice);
y=zeros(numThetan,nTau+1,nslice,NumElement);
Dy1=y;
Dy2=y;
Daligned=zeros(2*numThetan,numThetan*(nTau+1),nslice,NumElement);
for n=1:numThetan
    for ele=1:NumElement
        [y(n,:,:,ele),Dy1(n,:,:,ele),Dy2(n,:,:,ele)]=GaussianShift2D(project(:,:,ele,n),shift(:,n));
    end
    Daligned(2*n-1,n:numThetan:end,:,:)=Dy1(n,:,:,:);
    Daligned(2*n,n:numThetan:end,:,:)=Dy2(n,:,:,:);
end
gtot=zeros(N^2*NumElement,nslice);
f=0;
g_pert=zeros(1,2*numThetan);

for slice=1:nslice
    g=[];
    XTM=reshape(y(:,:,slice,:),numThetan*(nTau+1),NumElement);
    for ele=1:NumElement
        if(strcmp(frame,'EM'))
            thres=1;
            Mt=XTM(:,ele)+thres;
            Rdis=Ltol*xc(:,ele,slice)+thres;
            f=f+alpha(slice,ele)*sum(-log(Rdis).*Mt+Rdis);
            g= [g; alpha(slice,ele)*Ltol'*(-Mt./Rdis+1)];
            g_pert(1:2:end)=g_pert(1:2:end)-alpha(slice,ele)*log(Rdis)'*(Daligned(1:2:end,:,slice,ele))';
            g_pert(2:2:end)=g_pert(2:2:end)-alpha(slice,ele)*log(Rdis)'*(Daligned(2:2:end,:,slice,ele))';
        end
    end
    gtot(:,slice)=g;
end
g=[g_pert';gtot(:)];
