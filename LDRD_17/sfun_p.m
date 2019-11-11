function [f,g,XTM,f_sub]=sfun_p(shift,xc,XTM,Ltol) 
global frame N_delta  
global thetan numThetan dTau nTau 
global lambda sinoS N alpha
global nslice reg_str
%%======================Reconstruction with respect to shifts.
%%======================xc:N^2 x nchannel, current sample slices
nchannel=size(xc,2);
%%===========================================
[alignedSignal,Daligned]=GaussianShift1D(XTM,shift);
%%===========================================
XTM=reshape(alignedSignal,numThetan*(nTau+1),nchannel);

f=0;
g_pert=zeros(1,N_delta);
[~,~,Tik]=laplacian(N_delta);
for ele=1:nchannel
    if(strcmp(frame,'EM'))
        thres=1;
        Mt=XTM(:,ele)+thres;
        Rdis=Ltol*xc(:,ele)+thres;
        if(strcmp(reg_str,'TV'))
            epsilon=1e-4;
            f2_v=sqrt((Tik*shift).^2 + epsilon);
            f2=sum(sum(f2_v));
            g_reg= Tik'*Tik*shift./f2_v;
            % f2=1/2*sum((Tik*shift).^2);
            % g_reg=Tik'*Tik*shift;
            f=f+alpha(ele)*sum(-log(Rdis).*Mt+Rdis)+lambda*f2;
        end
        g_pert=g_pert-alpha(ele)*log(Rdis)'*(Daligned(:,:,ele))'+lambda*g_reg(:)';
    elseif(strcmp(frame,'LS'))
    end
end
g=g_pert';
