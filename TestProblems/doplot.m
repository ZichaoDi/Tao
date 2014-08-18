function doplot(v)
%--------------------------------------------------------------
% Plot results of soap-film minimal surface
%--------------------------------------------------------------
global y
n2 = length(v);
nm = sqrt(n2);
vs = reshape(v,nm,nm);
nx = nm + 2;
ny = nm + 2;
%--------------------------------------------------------------
% set up boundary conditions
%--------------------------------------------------------------
[x,y] = getborder(nx,ny);
%--------------------------------------------------------------
figure(3); 
%--------------------------------------------------------------
[xx,yy] = meshgrid(min(x):.02:max(x),min(y):.02:max(y));
[x,y] = meshgrid(x,y);
vv = interp2(x,y,vs,xx,yy,'cubic');
surfl(xx,yy,vv);
shading interp
colormap pink 
