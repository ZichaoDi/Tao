% Initialize paths
startup


% Note: This requires the Optimization Toolbox
%  GaussianFit1 needs fsolve for equation solving (only used for plotting)
% Note: This will clear the screen
XRFreconCompond

% OUTPUT: 
% XRF{i1,i2}(j1) = output from the i2 beamlet at the angle
%  thetan(i1) for j1st energy/detector channel (plot against DetChannel)
%
% Element{Z(j)} j=1:length(Z) shows element name
%


% So here is the data:
for i1=1:size(XRF,1)
    for i2=1:size(XRF,2)
        figure(100); 
      
semilogy(DetChannel,XRF{i1,i2}(:,1))
  title(sprintf('XRF at beam %d with angle %d',i2,thetan(i1)));  
  numphotons=sum(XRF{i1,i2}(:,1));
  xlabel('Energy Channel')
  xlim([min(DetChannel),max(DetChannel)])
  ylabel('Intensity')
  legend(['integral= ',num2str(numphotons)])
pause
    end
end



return

 UnitSpectra=rand(108,1000);
for i=1:108, UnitSpectra(i,:) = UnitSpectra(i,:)/sum(UnitSpectra(i,:)); end
plot(UnitSpectra(1,:))
hold on
plot(UnitSpectra(2,:),'g')

a=[];
for i=1:length(Z)
    PIXEL(i1,i2,:) = XRF{i1,i2}(i,2)*UnitSpectra(Z(i),:);
a = a + PIXEL(i1,i2,:);
end
% a tells me the cummulative spectra recorded by detector at i1,i2




% Stefan actually wants
%%OUT(i1,i2,1:1000) = spectra recorded for i1,i2





% % 
% % %%Input: 
% % Stefan wants to input:
% % % 
% % % for k = 1: E
% % %     W(i,j,k) = amount of element e in position i,j
% % % end
% %     
% % % Object and attentuation 
% % 
% % Output:


figure
clf;
for i=1:4
subplot(2,2,i)
imagesc(W(:,:,i));
colorbar
title(['Element ', num2str(i)])
end