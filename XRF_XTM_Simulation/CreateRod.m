%%================== Simulate glass rod sample for self-absorption 
xc=getCellCenteredGrid(omega,[N,N]);
xc=reshape(xc,N^2,2);
center = [1.3 -1.3 ; 0 0; 1.3 1.3].*Tol;
r = [ 0.05 1 0.05].*Tol;
if(NumElement==8)
    composition={'B','O','Na','Al','Si','K'};
    percentage=[0.26 2.111 0.07 0.04 0.81 0.01];
    Density=2.23;
    C=concentration(Density,percentage,composition);
else
    C=2.33;
end
w_ele=[19.3 C*120 19.3]; % units: g/cm^3
W = zeros(N^2,NumElement);
hollow=0;
for ele=1:NumElement
    if(ele==1)
        ind=1;
    elseif(ele==NumElement)
        ind=3;
    else
        ind=2;
    end
    pix = (xc(:,1)-center(ind,1)).^2+(xc(:,2)-center(ind,2)).^2 <= r(ind)^2;
    if(ele==1)
        W(find(pix==1),1)=w_ele(1);;
    elseif(ele==NumElement)
        W(find(pix==1),ele)=w_ele(ele);
    else
        if(hollow)
            pix = (xc(:,1)-center(2,1)).^2+(xc(:,2)-center(2,2)).^2 <= r(2)^2 & (xc(:,1)-center(2,1)).^2+(xc(:,2)-center(2,2)).^2 >= (2*r(2)/3)^2;
        end
        W(find(pix==1),ele)=w_ele(ele);
    end
end
W=reshape(W,[N,N,NumElement]);
% figure, for i=1:8,subplot(1,8,i);imagesc(W(:,:,i));colorbar;end
%%===================== Simulate on the reconstruction
load xstar_joint;
Wtemp=reshape(xstar_joint,N,N,3);
% Wtemp=reshape(x_admm(:,end),N,N,3);
smooth_rate=0;
W(:,:,1)=smooth2a(Wtemp(:,:,1),smooth_rate)/10;
W(:,:,6)=smooth2a(Wtemp(:,:,2),smooth_rate)/10;
W(:,:,8)=smooth2a(Wtemp(:,:,3),smooth_rate)/10;
% figure, for i=1:8,subplot(1,8,i);imagesc(W(:,:,i));colorbar;end
W=W/0.19;
%%=============================================

