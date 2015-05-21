%% Sven's sample
% subm=10;
% m=subm.*[3 3];
% NumElement=3;
W=zeros(m(1),m(2),NumElement);

class1x_1=1:1/3*m(1);
class1y_1=1:2/3*m(2);
class1x_2=1/3*m(1)+1:m(1);
class1y_2=1:1/3*m(2);
W(class1x_1,class1y_1,1)=0.1;
W(class1x_2,class1y_2,1)=0.1;

class2x_1=1/3*m(1)+1:2/3*m(1);
class2y_1=1/3*m(2)+1:2/3*m(2);
W(class2x_1,class2y_1,2)=0.5;


class3x_1=1:2/3*m(1);
class3y_1=2/3*m(2)+1:m(2);
class3x_2=2/3*m(2)+1:m(1);
class3y_2=1/3*m(2)+1:m(2);
W(class3x_1,class3y_1,3)=1;
W(class3x_2,class3y_2,3)=1;

% figure('name','sample');for i=1:NumElement, subplot(1,NumElement,i);imagesc(W(:,:,i));end
% figure, imagesc(sum(W,3));



