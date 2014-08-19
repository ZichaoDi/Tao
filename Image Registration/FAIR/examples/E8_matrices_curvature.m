%==============================================================================
% (c) Jan Modersitzki 2011/01/02, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% Tutorial for FAIR: regularization
%
% illustrates the tools for L2-norm based regularization
%
%  S(y) = alpha/2 * norm(B*(y-yRef)^2,
%
% where
%  alpha regularization parameter, weights regularization versus 
%        distance in the joint objective function, alpha = 1 here
%  yRef  is a reference configuration, e.g. yRef = x or a 
%        pre-registration result
%  B     a discrete partial differential operator either in explicit
%        matrix form or as a structure containing the necessary
%        parameters to compute B*y
% see also getCurvatureMatrix
%==============================================================================

clear, close all, help(mfilename);

                    % 2D example
omega  = [0,1,0,1]; % physical domin
m      = [16,12];   % number of discretization points
n      = prod(m);   % number of cells, and staggered dimensions
alpha  = 1;

% build curvature operator on a cell-centered and visualize
B = getCurvatureMatrix(omega,m);

% prepare the matrix free code
regularizer('reset','regularizer','mfCurvature','alpha',alpha);
[Sc, dS, d2S] = regularizer(zeros(2*n,1),omega,m);

% note: B*y can also be computed without storing B
%       the following loop computes A(:,j) + B(e(j)),
%       where mfBy is the matrix free operation B*y
%       and e(j) is the j-th unit vector;
%       hence A == B, (neither A nor B need to be stored)
e = @(j) ((1:size(B,2))' == j); % j-th unit vector
A = sparse(size(B,1),size(B,2));
for j=1:size(B,2), A(:,j) = d2S.By(e(j),omega,m); end;

% now the transpose
e = @(j) ((1:size(B,1))' == j); % j-th unit vector
C = sparse(size(B,2),size(B,1));
for j=1:size(B,1), C(:,j) = d2S.BTy(e(j),omega,m); end;

FAIRfigure(1); clf; 
subplot(2,3,1); spy(B);    title('B =curvature operator on cell-centered grid')
subplot(2,3,2); spy(A);    title('B from mfBy');
subplot(2,3,3); spy(B-A);  title('difference |B-mfBy()|');
subplot(2,3,4); spy(B');   title('B''')
subplot(2,3,5); spy(C);    title('B'' from mfBy');
subplot(2,3,6); spy(B'-C); title('difference |B''-mfBy(...,''BTy'')|');
