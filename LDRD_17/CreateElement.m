global cross_ind NumElement
if(strcmp(sample,'MRI'))
    A=double(imread('MRIhead-R.jpg'));
elseif(strcmp(sample,'Phantom'))
    A=phantom3d(N(1));
    cross_ind=floor(N(1)/2);
    A=squeeze(A(:,:,cross_ind));% sum(imread('phantom.png'),3);%
    A=map1D(double(A),[0,1]);
end
if(NumElement==1)
    W=sum(A,3);% abs(peaks(m(1)));%
else
    A(A(:)<0)=0;
    tol=eps^(1/2);
    W=zeros(m(1),m(2),NumElement);
        val=unique(A(:));
    if(length(val)>1)
        % val=val(val~=0);
        val_i=[];
        for i_v=1:length(val)
            val_i(i_v)=length(find(A(:)==val(i_v)));
        end
        ExtraVal=[0.1 0.2 0.3];
        [i1,i2]=sort(val_i,'descend');
        subind=[1 2 4];
        subind=subind(1:NumElement);
        if(length(val)<4)
            subind=[1 2];
            if(NumElement==1)
                subind=1;
            end
        end
        val=val(i2(subind));
        for i=1:length(subind)
            Ws=zeros(m(1),m(2));
            Ws(abs(A(:)-val(i))<tol)=val(i)+ExtraVal(i);
            W(:,:,i)=Ws;%+0.1;
        end
    end
end;

