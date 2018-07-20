function Output=gaussfilt2D(I,sigma);
M=size(I,1);
N=size(I,2);
[x,y]=meshgrid(1:M,1:N);
Output=zeros(size(I));
for i=1:M
    Output(i,:)=gaussfilt(1:N,I(i,:),sigma(1));
end
for i=1:N
    Output(:,i)=gaussfilt(1:M,Output(:,i),sigma(2));
end
% Exp_comp = -((x-mean([1:M])).^2+(y-mean([1:N])).^2)/(2*sigma(1)*sigma(2));
% Kernel= exp(Exp_comp)/(sqrt(2*pi)*sigma(1)*sigma(2));
% Output=conv2(I,Kernel','same');
