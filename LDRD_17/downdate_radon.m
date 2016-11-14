function [LH,nThe_H,nT_H]=downdate_radon(L,res_prob,numThetan,nTau)
%%===== Directly construct coarsened version of radon transform L, rows of L correspond to measurements, column of L corresponds to variables. Rows are downdated by taking every even number of beamlet and then take a full weighting on the angle. Columns are downdated by full weighting.
global N level 
sub_tau=2:2:nTau+1;
mtol=size(L,2);
if(res_prob==0)
    L=reshape(full(L),numThetan,(nTau+1),mtol);
    for i=1:numThetan
        L_resTau(i,:,:) = restrict_residual(nTau+1)*squeeze(L(i,:,:));
    end
    L=L_resTau;
    nT_H=size(L,2);
    % if(numThetan>1)
    %     LH_left=restrict_residual(numThetan)*reshape(L,numThetan,nT_H*mtol);
    % else
        LH_left=reshape(L,numThetan,nT_H*mtol);
    % end

    nThe_H=size(LH_left,1);
    LH_left=reshape(LH_left,nThe_H*nT_H,mtol);
    LH=zeros(nThe_H*nT_H,N(level)^2);
    for i=1:nThe_H*nT_H
        LH(i,:)=downdate(LH_left(i,:),1);
    end

elseif(res_prob==1)
    L=reshape(L,numThetan,nTau+1);
    for i=1:numThetan
        L_resTau(i,:) =restrict_residual(nTau+1)*L(i,:)';
    end
    % if(numThetan>1)
    %     LH_left=restrict_residual(numThetan)*L_resTau;%L(:,sub_tau);
    % else
        LH_left=L_resTau;%L(:,sub_tau);
    % end
    nThe_H=size(LH_left,1);
    nT_H=size(LH_left,2);
    LH=reshape(LH_left,nThe_H*nT_H,1);
end
nT_H=nT_H-1;
