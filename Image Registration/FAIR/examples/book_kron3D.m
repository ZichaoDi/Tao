%$ (c) Jan Modersitzki 2009/04/07, see FAIR.2 and FAIRcopyright.m.
%$ \url{http://www.cas.mcmaster.ca/~fair/index.shtml}

% memory efficient implementation of
% reshape(Mx,m), Mx = kron(M{1},kron{M{2},M{3})) * x(:)

x = reshape(x,m); % make sure x is a 3d-array of size m

% given a permutation L of [1,2,3], the following function
% 1. permutes x, such that j=J(1) is the first dimension
% 2. reshapes the permuted x, such that is m(j)-by-M(J(2))*m(J(3))
% 3. multiplies this by M which is assumed to be m(j)-by-m(j)
% 4. undoes the reshape to make the result m(j)-by-m(J(2))-ny-m(J(3))
% 5. undoes the permute

operate = @(M,x,L) ipermute( reshape( ... 
  M*reshape( permute(z,L), m(L(1)),[]), m(L)), L);
% run over all directions
for ell=1:3,
  % make the ell-th component the first 
  L = [ell,setdiff(1:3,ell)];  
  % operate as indicated above
  x = operate(M{ell},x,L);
end;