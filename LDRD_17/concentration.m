function C = concentration(density,percentage,composition);
load PeriodicTable
NA=6.02e23;%Avogadro's number
if(ismac)
     loadlibrary('/opt/local/lib/libxrl.dylib','/opt/local/include/xraylib/xraylib.h');
else
     loadlibrary('/homes/wendydi/local/xraylib/lib/libxrl.so','/homes/wendydi/local/xraylib/include/xraylib/xraylib.h');
end
%%======================================================================
subElement=length(composition);
for i=1:subElement
    Z(i)=find(strcmp(Element,composition{i}));
    ElementDensity(i)=calllib('libxrl','ElementDensity',Z(i));%
    A(i)=calllib('libxrl','AtomicWeight',Z(i));
end
A_bar=sum(percentage.*A);
na=density*NA/A_bar.*percentage;
C=na.*A/NA;

