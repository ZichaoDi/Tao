
function [specP]=GaussianFit1(Energy,Intensity)
global DetChannel numChannel plotSpecSingle NumElement

p=zeros(NumElement,1,numChannel);
specP=zeros(numChannel,1);
if(plotSpecSingle)
    figure('name','XRF elemental spectrum')
end
for i=1:NumElement
    
    if(Intensity(i)>1e-10)
        miu=Energy(i);
        theta=Intensity(i)^(1/3);
        % fprintf('%i: mean = %d, standard deviation = %d \n', i, miu,theta);
        p(i,1,:)=Intensity(i)/(theta*sqrt(2*pi))*exp(-(DetChannel-miu).^2./(2*theta^2));
        if(plotSpecSingle)
            plot(DetChannel,reshape(p(i,1,:),numChannel,1),'r-')
            title('Gaussian spectrum for each element')
            hold on;
        end
        specP=specP+reshape(p(i,1,:),numChannel,1);
    end
end
