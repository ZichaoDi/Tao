global numThetan maxiter W0 sinoS
global totN xiter fiter ErrIter alpha

ReconAttenu = 1; % 0: Recover W; 1: Recover miu
slice=1;
do_setup;
rng('default');
res=3;
d0=[0 -res -res 0     res res res 0   -res;...
    0  0   -res -res -res 0   res res res];
d0=d0(:,3);
initial_direction=size(d0,2);

maxiter=200;
recon=sparse(N^2*NumElement*nslice+numThetan*2,1);
totN=N^2*NumElement*nslice+numThetan*2;
W0=recon;
Mt(isnan(Mt))=0;
Mt(isinf(Mt))=0;
Mt=Mt(:,:,1,:);
alpha=zeros(nslice,NumElement);
for slice=1:nslice
    for ele=1:NumElement
        alpha(slice,ele)=std(reshape(Mt(:,slice,ele,:),(nTau+1)*numThetan,1));
    end
    alpha(slice,:)=alpha(slice,:)./sum(alpha(slice,:));
end

    NF = [0*N; 0*N; 0*N];
    Lmap=sparse(L);
    for res_step=1:initial_direction
        x0 = zeros(totN,1);
        err0=norm(W0-x0);
        fctn_COR=@(x)sfun_shift_xy(x,Mt,Lmap);% on attenuation coefficients miu;
        low=zeros(totN,1);
        up=inf*ones(totN,1);
        low(1:2:numThetan*2)=-nTau/2*ones(numThetan,1);
        low(2:2:numThetan*2)=-nslice/2*ones(numThetan,1);
        up(1:2:numThetan*2)=nTau/2*ones(numThetan,1);
        up(2:2:numThetan*2)=nslice/2*ones(numThetan,1);
        bounds=1;
        if(bounds)
            [xCOR,f,g,ierror] = tnbc (x0,fctn_COR,low,up); % algo='TNbc';
        else
            [xCOR,f,g,ierror] = tn (x0,fctn_COR);
        end
    end
