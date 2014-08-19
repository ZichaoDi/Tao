function  [T,R]=CreateSimpleTests(c,m,r);

R=zeros(m(1),m(2));
T=zeros(m(1),m(2));
for i=1:m(1)
    for j=1:m(2)
        if((i-c(1))^2+(j-c(1))^2<=r^2)
            T(i,j)=255;
        else
            T(i,j)=1;
        end
        if((i-c(2))^2+(j-c(2))^2<=r^2)
                R(i,j)=255;
        else
            R(i,j)=1;
        end
    end
end

% close all
% imshow(T)
% figure,imshow(R)
% axis xy
