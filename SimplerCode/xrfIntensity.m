% Intensity=IncidentFlux*IlluminationArea*CrossSection*FluorescenceYield*(1-SelfAbsorption)
function [xrf_intensity,beam_energy]=xrfIntensity
close all
clc
warning off;
options=optimset('Display','off');
loadlibrary('/opt/local/lib/libxrl.dylib','/opt/local/include/xraylib/xraylib.h');
% Z=[30 33 34 20 22 23  26 27 28 29 35 24 25]; %% elements
Z=[30 20 35 26];
Line=2; %% Transition Line
flux=1e9;
G=1e-5;
xrf_intensity=[];
FWHM=0;
i=1;
beam_energy=20;
new_energy=beam_energy;
n=100;
fprintf(1,'Element     Line     FluorescenceEnergy     Intensity\n')
while(i<=length(Z))
    density=calllib('libxrl','ElementDensity',Z(i));
    thickness=0.1;
    mu_0=calllib('libxrl','CS_Total',Z(i),beam_energy);
    j=1;
    while (j<=length(Line))
        new_energy=calllib('libxrl','LineEnergy',Z(i),-Line(j));
        mu_1=calllib('libxrl','CS_Total',Z(i),new_energy);
        Q=calllib('libxrl','CS_FluorLine_Kissel',Z(i),-Line(j),beam_energy);
        w=calllib('libxrl','FluorYield',Z(i),Line(j)-1);
        chi=mu_0+mu_1;
        A_corr=1-exp(-chi*density*thickness);
        intensity=flux*G*Q*w*A_corr;
        old_energy=new_energy;
        
        if(intensity~=0)
            fprintf('    %i        %i          %5.2f             %5.2f\n', Z(i), -Line(j), new_energy, intensity);
        x0=new_energy;
        theta=20/intensity^(1/4);
        [miu]=fsolve(@(miu)GaussianFit(miu,theta,intensity,FWHM,n),x0,options);
        %%%%%======= visualize the spectrum obeying a normal distribution
        ind1=miu+FWHM;
        FWHM=FWHM+10*2.354*beam_energy;
%         x=linspace(ind1-4*2.354*theta,ind1+4*2.354*theta,n);
%         p=1/(theta*sqrt(2*pi))*exp(-(x+theta-ind1-theta).^2./(2*theta^2));
%         plot(x,p,'r-')
%         hold on;
        else
            miu=0;
            theta=0;
        end
        xrf_intensity(length(Line)*(i-1)+j,1:4)=[new_energy,intensity,miu,theta];

        %%%%===================================================================
        j=j+1;
    end
    fprintf('=====================================================\n')
    i=i+1;
end


function err=GaussianFit(miu,theta,I,FWHM,n)
ind1=miu+FWHM;
t=linspace(ind1-4*2.354*theta,ind1+4*2.354*theta,n);
p=1/(theta*sqrt(2*pi))*exp(-(t+theta-ind1-theta).^2./(2*theta^2));
I0=sum(p);
err=norm(I-I0);
        







% [temp,upind]=sort(xrf_intensity(:,1));
% xrf_intensity=xrf_intensity(upind,:);
% for i=1:size(xrf_intensity,1)
%     if(i==1)
%         psi1=(xrf_intensity(i,1)-0)/2;
%         psi2=(xrf_intensity(i+1,1)-xrf_intensity(i,1))/2;
%         psi=min(psi1,psi2);
%     elseif(i==size(xrf_intensity,1))
%         psi1=(xrf_intensity(i,1)-xrf_intensity(i-1,1))/2;
%         psi2=(beam_energy-xrf_intensity(i,1))/2;
%         psi=min(psi1,psi2);
%     else
%         psi1=(xrf_intensity(i,1)-xrf_intensity(i-1,1))/2;
%         psi2=(xrf_intensity(i+1,1)-xrf_intensity(i,1))/2;
%         psi=min(psi1,psi2);
%     end
%
%     FWHM=1/2*2.354*psi; %psi: width;
%     x=linspace(xrf_intensity(i,1)-FWHM,xrf_intensity(i,1)+FWHM,100);
%     p=xrf_intensity(i,2)/sqrt(2*pi*psi^2)*exp(-(x-xrf_intensity(i,1)).^2./(2*psi^2));
%     plot(x,p,'r-')
%     hold on;
% end
% % figure,plot(xrf_intensity(:,1),xrf_intensity(:,2),'r.-')



