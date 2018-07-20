global numThetan maxiter W0 sinoS
global totN xiter fiter ErrIter alpha

ReconAttenu = 0; % 0: Recover W; 1: Recover miu
slice=1;
do_setup;
rng('default');
res=3;
d0=[0 -res -res 0     res res res 0   -res;...
    0  0   -res -res -res 0   res res res];
d0=d0(:,3);
initial_direction=size(d0,2);

maxiter=400;
totN=N^2*NumElement*nslice+numThetan*2;
W0=sparse(totN,1);
for n=1:numThetan
    dd=Mt(:,:,:,n);
    dd(isnan(dd))=min(dd(:));
    Mt(:,:,:,n)=dd;
end 
Mt(isinf(Mt))=0;
Mt=Mt(:,1:nslice,1:NumElement,:);
alpha=ones(nslice,NumElement);
for slice=1:nslice
    for ele=1:NumElement
        alpha(slice,ele)=std(reshape(Mt(:,slice,ele,:),(nTau+1)*numThetan,1));
        % temp=Mt(:,slice,ele,:);
        % Mt(:,slice,ele,:)=temp./max(temp(:))+abs(min(temp(:)));
    end
    alpha(slice,:)=alpha(slice,:)./sum(alpha(slice,:));
end
Mt=Mt./max(Mt(:))+abs(min(Mt(:)));%%==normalize data;
alpha=ones(size(alpha));

NF = [0*N; 0*N; 0*N];
Lmap=sparse(L);
MtPad=GaussianPadded(Mt);
for res_step=1:initial_direction
    load opt3D_xy_normalized;
    x0 = xCOR;%zeros(totN,1);
    err0=norm(W0-x0);
    fctn_COR=@(x)sfun_shift_xy(x,MtPad,Lmap);% on attenuation coefficients miu;
    low=zeros(totN,1);
    up=inf*ones(totN,1);
    low(1:2:numThetan*2)=floor(-nTau/2*ones(numThetan,1));
    up(1:2:numThetan*2)=floor(nTau/2*ones(numThetan,1));
    low(2:2:numThetan*2)=floor(-nslice/2*ones(numThetan,1));
    up(2:2:numThetan*2)=floor(nslice/2*ones(numThetan,1));
    bounds=1;
    x0(1:2:n_delta)=-align(ang_ind,1);
    x0(2:2:n_delta)=-align(ang_ind,2);
    [~,~,aligned]=feval(fctn_COR,x0);
    fctn=@(x)sfun3D(x,aligned,Lmap);% on attenuation coefficients miu;
    W0=W0(n_delta+1:end);
    [x,f,g]=tnbc(x0(n_delta+1:end),fctn,low(n_delta+1:end),up(n_delta+1:end));
    % return;
    if(bounds)
        [xCOR,f,g,ierror] = tnbc (x0,fctn_COR,low,up); % algo='TNbc';
    else
        [xCOR,f,g,ierror] = tn (x0,fctn_COR);
    end
end

save('opt3D_xy_normalized.mat','xCOR');
