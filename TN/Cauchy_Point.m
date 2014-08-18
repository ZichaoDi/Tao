function [x_cp,i_cauchy,alpha_cp,p_cp]=Cauchy_Point(x0,g, accrcy, xnorm, sfun,low,up)

x_cp=x0;
i_cauchy=0;
alphaDis=zeros(size(x0));
indl=find(g>0);
indu=find(g<0);
alphaDis(indl)=(x0(indl)-low(indl))./g(indl);
alphaDis(indu)=(x0(indu)-up(indu))./g(indu);
alphaDis=min(alphaDis,1);
alpha=[0;sort(unique(alphaDis(alphaDis~=0)))];
p_cp=[];
alpha_cp=[];
for i=1:length(alpha)-1
    p=-g;
    p(alpha(i+1)>alphaDis)=0;
    p_cp=[p_cp,p];
    [~, ~, x] = crash (x0-alpha(i)*g, low, up);
    gvp = gtims (p, x0, g, accrcy, xnorm, sfun);
    gt=g'*p+x'*gvp;
    Ht=p'*gvp;
    if(Ht<0)
        disp('Hessian is negative')
    end
    delta=-gt/Ht;
    alpha_cp=[alpha_cp,alpha(i)];
    if(gt>0 & i~=1)
        p_cp=p_cp(:,1:end-1);
        x_cp=x;
        i_cauchy=1;
        break;
    elseif(i<length(alpha))

        if(delta>=0 & delta<alpha(i+1)-alpha(i))
            x_cp=x+delta*p;
            i_cauchy=-1;
            alpha_cp=[alpha_cp,delta];
            break;
        end
    end
    
end
