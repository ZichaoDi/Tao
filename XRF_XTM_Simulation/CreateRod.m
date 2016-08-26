%%================== Simulate glass rod sample for self-absorption 
xc=getCellCenteredGrid(omega,[N,N]);
xc=reshape(xc,N^2,2);
% center = [1.4 -1.2 ; 0 0; 1.1 1.1].*Tol;
center = [2.9 -2.7 ; 0 0; 2.6 2.6].*Tol;
r = [ 0.17 3.0 0.17].*Tol;
w_ele=[19.3 2.33 19.3]; % units: g/cm^3
w_ele=w_ele.*[1 0.5 1]*2;
W = zeros(N^2,3);
hollow=0;
for ele=1:NumElement
    pix = (xc(:,1)-center(ele,1)).^2+(xc(:,2)-center(ele,2)).^2 <= r(ele)^2;
    if(ele==1)
        % W(find(pix==1),3)=w_ele(3);%max(reshape(iR(:,:,3),N^2,1));
        W(find(pix==1),1)=w_ele(1);%max(reshape(iR(:,:,1),N^2,1));
        % W(find(pix==1),2)=w_ele(2)+1;%6e-3;
    elseif(ele==2)
        % rec=[[1.1e-5 1.1e-5+r(3)];[1.2e-5 -1.2e-5+r(3)];[1.1e-5 1.1e-5-r(3)];[1.2e-5 -1.2e-5-r(3)]];
        % rec=[center(1,:);center(1,:)-[r(3)/2,0];center(3,:)-[r(3)/2,0];center(3,:)];
        % pix_extra=isPointInPolygon(xc,rec);
        if(hollow)
            pix = (xc(:,1)-center(ele,1)).^2+(xc(:,2)-center(ele,2)).^2 <= r(ele)^2 & (xc(:,1)-center(ele,1)).^2+(xc(:,2)-center(ele,2)).^2 >= (2*r(ele)/3)^2;
        end
        W(find(pix==1),2)=w_ele(2);%6e-3;
        % W(find(pix_extra==1),2)=w_ele(2);%6e-3;
    elseif(ele==NumElement)
        W(find(pix==1),ele)=w_ele(ele);%max(reshape(iR(:,:,ele),N^2,1));
        % W(find(pix==1),2)=w_ele(2)+1;%6e-3;%max(reshape(iR(:,:,2),N^2,1));
    end
end
W=reshape(W,[N,N,3]);
