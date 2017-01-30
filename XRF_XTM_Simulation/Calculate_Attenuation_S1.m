function [ConstSub, f]=Calculate_Attenuation_S1(W,NumElement,L,GlobalInd,SelfInd,m,nTau,xrfData,Mt, M)
%%==== Given elemental map W and pre-calculate the beam and fluorescent attenuation coefficients
global numChannel NoSelfAbsorption numThetan
global MU_e area_xrf
global TempBeta Beta Joint frame
mtol=prod(m);
W=reshape(W,mtol,NumElement);
MU=zeros(prod(m),NumElement);
for i=1:NumElement
    temp=sum(reshape(W,prod(m),NumElement)*MU_e(:,:,i+1),2);
    temp=flipud(reshape(temp,m(1),m(2))');
    MU(:,i)=temp(:);
end
%%%%% ====================================================================
if(norm(W(:))>0)
InTens=ones(numThetan*(nTau+1),mtol);
OutTens=ones(numThetan*(nTau+1),mtol,NumElement);
    for n=1:numThetan
        for i=1:nTau+1
            ind_bt=(i-1)*numThetan+n;
            index=GlobalInd{ind_bt};
            if(~isempty(index))
                index_sub=sub2ind(m,index(:,2),index(:,1));
                msub=length(index_sub);
                for v_count=1:msub
                    v=index_sub(v_count);
                    if(~isempty(SelfInd{ind_bt,v}{1}))
                        InTens(ind_bt,v)=exp(-sum(sum(W(SelfInd{ind_bt,v}{1},:).*SelfInd{ind_bt,v}{3})));
                    end
                    if( ~isempty(SelfInd{ind_bt,v}{2}) && ~NoSelfAbsorption)
                        OutTens(ind_bt,v,:)=exp(-sum(MU(SelfInd{ind_bt,v}{2},1:NumElement),1)./(length(SelfInd{ind_bt,v}{2})+1)*area_xrf(ind_bt,v));
                    end
                end
            end
        end
    end
% save('testC.mat','L','InTens','OutTens','M');
tic;
ConstSub=bsxfun(@times,reshape(bsxfun(@times,full(L.*InTens),OutTens),[numThetan*(nTau+1),1,mtol,NumElement]),reshape(M',[1,numChannel,1, NumElement])); 
toc;
size(ConstSub)
% save ConstSub ConstSub
% L=full(L);
% I=full(InTens);
% O=OutTens;
% fidL=fopen('L.txt','a');
% fidI=fopen('I.txt','a');
% fidO=fopen('O.txt','a');
% [N1,N2]=size(L);
% [N1,N2,N3]=size(O);
% fprintf(fidL,'%d    %d     \n',N1,N2);
% fprintf(fidI,'%d    %d     \n',N1,N2);
% fprintf(fidO,'%d    %d    %d \n',N1,N2,N3);
% for i=1:N1
%     for k=1:N2
%         if(L(i,k)~=0)
%             fprintf(fidL,'%d    %d     %f\n',i,k,L(i,k)');
%         end
%         if(I(i,k)~=0)
%             fprintf(fidI,'%d    %d     %f\n',i,k,I(i,k)');
%         end
%         for l=1:N3
%             fprintf(fidO,'%d    %d     %d   %f\n',i,k,l,O(i,k,l)');
%         end
%     end
% end
% 
% 
% fidM=fopen('M.txt','a');
% [N4,N5]=size(M);
% fprintf(fidM,'%d    %d     \n',N4,N5);
% for l=1:N4
%     for j=1:N5
%         if(M(l,j)~=0)
%             fprintf(fidM,'%d    %d     %f\n',l,j,M(l,j)');
%         end
%     end
% end
% 
% fclose(fidI),fclose(fidL);fclose(fidO);fclose(fidM);
clear M InTens OutTens
return
else
    ConstSub=bsxfun(@times,reshape(full(L),numThetan*(nTau+1),1,mtol,1),reshape(M',[1,numChannel,1,NumElement]));

end
ConstSub=sparse(reshape(ConstSub,[numThetan*(nTau+1)*numChannel,mtol*NumElement]));
% ConstSub=kron(M',L.*InTens.*squeeze(sum(OutTens,3)./NumElement)); %% Igonor self-absorption
% clear OutTens InTens
% ConstSub=kron(ones(numChannel,1),kron(ones(1,NumElement),reshape(L.*InTens,numThetan*(nTau+1),mtol)).*reshape(OutTens,numThetan*(nTau+1),mtol*NumElement)).*kron(M',ones(numThetan*(nTau+1),mtol)); %% Kronecker version 
%#############################################
%%================== Memory saving version
% c_1=reshape(bsxfun(@times,full(L.*InTens),OutTens),[numThetan*(nTau+1),mtol*NumElement]);
% ConstSub=sparse(numThetan*(nTau+1)*numChannel,mtol*NumElement);
% for i1=1:numChannel
%     for i2=1:numThetan*(nTau+1)
%           ConstSub(i2+(i1-1)*numThetan*(nTau+1),:)=c_1(i2,:).*M(i1,:); % reshape(bsxfun(@times,squeeze(c_1(i2,:,:)),M(:,i1)'),mtol*NumElement,1);% 
%     end
% end
% clear M InTens OutTens
%%========================================================
XRF_v=ConstSub*W(:);
if(Joint==0)
    if(strcmp(frame,'EM'))
        thres=1;
        f=sum(-log(XRF_v+thres).*(xrfData(:)+thres)+XRF_v+thres);
    else
        f=sum((XRF_v-xrfData).^2);  
    end
elseif(Joint==1)
    MU_XTM=W*squeeze(MU_e(:,1,1));
    Rdis=reshape(L,numThetan*(nTau+1),mtol)*MU_XTM;
    if(strcmp(frame,'EM'))
        thres=1;
        f=sum(-log(XRF_v+thres).*(xrfData(:)+thres)+XRF_v+thres);
        f_XTM=sum(-log(Rdis+thres).*(Mt+thres)+Rdis+thres);
    elseif(strcmp(frame,'LS'))
        f=sum((XRF_v-xrfData(:)).^2);  
        f_XTM=sum((Rdis-Mt).^2);
    elseif(strcmp(frame,'mix'))
        thres=1;
        f=sum(-log(XRF_v+thres).*(xrfData(:)+thres)+XRF_v+thres);
        f_XTM=sum((Rdis-Mt).^2);
    end
    f = TempBeta*f + Beta*f_XTM;
end
