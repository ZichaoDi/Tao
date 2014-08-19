%$ function y = rigid2D(w,x)
%$ (c) Jan Modersitzki 2009/04/07, see FAIR.2 and FAIRcopyright.m.
%$ \url{http://www.cas.mcmaster.ca/~fair/index.shtml}
function yc = rigid2Dbook(wc,xc)

xc = reshape(xc,[],2);
yc = [(cos(wc(1))*xc(:,1) - sin(wc(1))*xc(:,2) + wc(2));
      (sin(wc(1))*xc(:,1) + cos(wc(1))*xc(:,2) + wc(3))];
