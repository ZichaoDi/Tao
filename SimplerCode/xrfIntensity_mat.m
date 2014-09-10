% Intensity=IncidentFlux*IlluminationArea*CrossSection*FluorescenceYield*(1-SelfAbsorption)
%function [xrf_intensity,beam_energy]=xrfIntensity_mat
global DetChannel numChannel
close all
clc
warning off;

% Z=[30 33 34 20 22 23  26 27 28 29 35 24 25 50 60 70 80]; %% elements
Z=[30 20 35 26];
load PeriodicTable
% disp('Existing elements:'),disp({Element{Z}});
Line=2; %% Transition Line
load ElementStrKalpha2
flux=1e5;
G=1e-4;
xrf_intensity=[];
FWHM=0;

i=1;
beam_energy=20;
new_energy=beam_energy;
numChannel=1e3;
DetChannel=linspace(0,100,numChannel)';
netSpec=zeros(numChannel,1);
fprintf(1,'Element     AtomicNumber    Line     FluorescenceEnergy     Intensity\n')
while(i<=length(Z))
    Density=density(Z(i));
    thickness=0.01;
    mu_0=mu0(Z(i));
    j=1;
    while (j<=length(Line))
        new_energy=Kalpha2Energy(Z(i));
        mu_1=mu1(Z(i));
        Q=CS(Z(i));
        w=1;%FluorYield(Z(i));
        chi=mu_0+mu_1;
        A_corr=1-exp(-chi*Density*thickness);
        intensity=flux*G*Q*w*A_corr;
        old_energy=new_energy;
        
        if(intensity~=0)
            fprintf('    %s           %d          %i          %5.2f               %5.2f\n', Element{Z(i)}, Z(i),-Line(j), new_energy, intensity);
        miu=new_energy;
        theta=Intensity(i)^(1/3);%1e4/Intensity(i)/sqrt(2*pi);
%          fprintf('%i: mean = %d, standard deviation = %d \n', i, miu,theta);
        p(i,1,:)=Intensity(i)/(theta*sqrt(2*pi))*exp(-(DetChannel-miu).^2./(2*theta^2));
        else
            miu=new_energy;
            theta=0;
        end
        xrf_intensity(length(Line)*(i-1)+j,1:4)=[new_energy,intensity,miu,theta];
        
        %%%%===================================================================
        j=j+1;
    end
    netSpec=netSpec+p;
    fprintf('=====================================================\n')
    i=i+1;
end
I(1,1)=mu0(Z(1));
I(1,2)=mu0(Z(2));
I(1,3)=mu0(Z(1));
I(2,1)=mu0(Z(2));
I(2,2)=0.5*(mu0(Z(3))+mu0(Z(4)));
I(2,3)=mu0(Z(2));
I(3,1)=mu0(Z(3));
I(3,2)=mu0(Z(1));
I(3,3)=mu0(Z(2));
I=I.*1e-1;
% plot(DetChannel,netSpec,'r-')

% function err=GaussianFit(miu,theta,I,FWHM,n)
% ind1=miu+FWHM;
% t=linspace(ind1-4*2.354*theta,ind1+4*2.354*theta,n);
% p=1/(theta*sqrt(2*pi))*exp(-(t+theta-ind1-theta).^2./(2*theta^2));
% I0=sum(p);
% err=norm(I-I0);