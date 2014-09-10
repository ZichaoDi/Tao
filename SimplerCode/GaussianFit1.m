
function [specP]=GaussianFit1(Energy,Intensity)
global DetChannel numChannel
num=length(Energy);
FWHM=100;

p=zeros(num,1,numChannel);
specP=zeros(numChannel,1);
% figure('name','Individual Elemental Spectrum')
for i=1:num
    
    if(Intensity(i)>1e-10)
        %
        % theta=20/Intensity(i)^(1/4);
        % options=optimset('Display','off');
        % [miu]=fsolve(@(miu)GaussianFit(miu,theta,Intensity(i),FWHM,n),Energy(i),options);
        % %%%%%======= visualize the spectrum obeying a normal distribution
        % ind1=miu+FWHM;
        %
        % FWHM=FWHM+2.354*10*beam_energy;
        % x=linspace(ind1-4*2.354*theta,ind1+4*2.354*theta,numChannel);
        % p=1/(theta*sqrt(2*pi))*exp(-(x+theta-ind1-theta).^2./(2*theta^2));
        % feas=find(x>0);
        % plot(x(feas),p(feas),'r.-');hold on;
        miu=Energy(i);
        theta=Intensity(i)^(1/3);%1e4/Intensity(i)/sqrt(2*pi);
%          fprintf('%i: mean = %d, standard deviation = %d \n', i, miu,theta);
        p(i,1,:)=Intensity(i)/(theta*sqrt(2*pi))*exp(-(DetChannel-miu).^2./(2*theta^2));
        
%         plot(DetChannel,reshape(p(i,1,:),numChannel,1),'r-')
%         hold on;
        specP=specP+reshape(p(i,1,:),numChannel,1);
    end
end




function err=GaussianFit(miu,theta,I,FWHM,n)
ind1=miu+FWHM;
t=linspace(ind1-4*2.354*theta,ind1+4*2.354*theta,n);
p=1/(theta*sqrt(2*pi))*exp(-(t+theta-ind1-theta).^2./(2*theta^2));
I0=sum(p);
err=norm(I-I0);
