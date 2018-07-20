function y = dirac(s,n)

% dirac - dirac function.
%
%   y = dirac(s,n);
%
%   s is the size of the matrix, and 
%   n is the location of the dirac (default (1,...,1)).
%
%   Copyright (c) 2004 Gabriel Peyré

if nargin<2
    n = ones(size(s));
end

if length(s)==1
    y = zeros(s,1);
else
    y = zeros(s);
end
y = array_set_val( y, n, 1 );

function y = array_set_val( x, ind, v )

% array_set_val - set a value in a multidimensional array using a vector for the index
%
%   y = array_set_val( x, ind, v );
%
%   Copyright (c) 2004 Gabriel Peyré

ind = num2cell(ind);
x( ind{:} ) = v;
y = x;
