omega = [0 1 0 1 0 1]; m = [5 4 6];
dim   = 3
omega = omega(1:2*dim), m = m(1:dim)
yC    = getCellCenteredGrid(omega,m);
e     = @(i) (1:length(yC)==i)';

fctn = @stg2center;

P = fctn(m);
A = zeros(fliplr(size(P)));

for i=1:size(A,2),
  A(:,i) = fctn(e(i),m);
end;

figure(1); clf; spy(A-P')

