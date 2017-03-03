function [LH]=downdate_radon(L,numThetan,nTau)
%%===== Directly construct coarsened version of radon transform L, rows of L correspond to measurements, column of L corresponds to variables. Rows are downdated by taking every even number of beamlet and then take a full weighting on the angle. Columns are downdated by full weighting.
global N_level  
LH=zeros(numThetan*(nTau+1),N_level(2)^2);
for i=1:numThetan*(nTau+1)
    LH(i,:)=downdate(L(i,:),1);
end

