%$ function yc = affine2D(wc,xc)
%$ (c) Jan Modersitzki 2009/04/06, see FAIR.2 and FAIRcopyright.m.
%$ \url{http://www.cas.mcmaster.ca/~fair/index.shtml}
function yc = affine2DBook(wc,xc)

xc = reshape(xc,[],2);
yc = [(wc(1)*xc(:,1) + wc(2)*xc(:,2) + wc(3))
      (wc(4)*xc(:,1) + wc(5)*xc(:,2) + wc(6))];
