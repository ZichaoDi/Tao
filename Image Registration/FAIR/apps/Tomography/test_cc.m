clear all
close all


load CoarsedT
load CoarsedR
T=1./CoarsedT;
R=1./CoarsedR;
inter('reset','inter','linearInter','regularizer','moments','theta',1e-1);
omega=[0 size(R,1) 0 size(R,2)];
m = floor(size(R)/1);
T = inter('coefficients',T,[],omega,'regularizer','moments','theta',1e1);
R = inter('coefficients',R,[],omega,'regularizer','moments','theta',1e1);
T=-log(T);
R=-log(R);
figure, viewImage2D(GrayScale(T),omega,m);
figure, viewImage2D(GrayScale(R),omega,m);

 Q = ifft2(fft2(R-mean(R(:))) .* conj(fft2(T-mean(T(:)))));

Q=sqrt(real(Q).^2+imag(Q).^2);
[v1, I2] = max(flipud(Q));
[v2,I1]=max(v1);
I2=I2(I1);
peak=[I1,I2]
Q(peak(1),peak(2))
figure,
imagesc(Q)
colormap jet