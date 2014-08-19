function Tc=GrayScale(T)

Tc=double(T(:));
minTc=min(Tc);maxTc=max(Tc);dG=maxTc-minTc;
dG=dG+2*(dG==0);
Tc=255/double(dG)*double(Tc-minTc);
Tc=reshape(Tc,size(T));