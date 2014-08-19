% (c) Jan Modersitzki 2010/12/27, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% function Qw = matVecQw(w,Q,flag)
%
% computes y = Q*w or Q'*w
% where Q =kron(Q3,kron(Q2,Q1)) is efficiently stored as Q= {Q1,Q2,Q3}
% see also splineTransformation2Dsparse 

function Qw = matVecQw(w,Q,flag)

if nargin == 0,
  help(mfilename);
  return;
end;

% default is Q*w
if ~exist('flag','var'), flag = 'Qw';  end;

if strcmp(flag,'QTw'),
  % build transposes, note: small matrices
  for i=1:length(Q),    Q{i}=Q{i}';  end;
  Qw = matVecQw(w,Q);
  return;
end;

% reconstruct sizes m and p from Q and w, Q{j} is m{j}-by-p{j}

dim = length(Q);
for j=1:length(Q),
  p(j) = size(Q{j},2);
  m(j) = size(Q{j},1);
end;

switch dim,
  case 1,
    Qw = Q{1}*w;
  case 2,
    w  = reshape(w,[p,dim]);
    Qw = zeros(m(1),m(2),dim);
    for i =1:dim,
      Qw(:,:,i) = Q{1}*w(:,:,i)*Q{2}'; 
    end;
  case 3,
    w  = reshape(w,[p,dim]);
    Qw = zeros(m(1),m(2),m(3),dim);
    for i =1:dim,
        Qw(:,:,:,i) = tensorProd(Q,w(:,:,:,i));
    end;    
end;

Qw = reshape(Qw,[],1);

return;
% reshape w to fit p
w = reshape(w,p);

if length(p) == 2, 
  % the 2D case: kron(Q2,Q1)*w(:) = (Q1*w*Q2')(:)
  w = Q{1}*w*Q{2}';
else
  % the 3D case, note, dim-1 needs no reorganization and is thus handled at last
  % 3rd DIR: make 3-dim leading, reshape w to p(3)-by-p(1)*p(2),
  % multiply by Q{3} and reorganize result
  w = permute(w, [3 1 2]);
  w = reshape(w,p(3),p(1)*p(2));
  w = Q{3}*w;
  w = reshape(w,m(3),p(1),p(2));
  w = permute(w,[2 3 1]);

  % 2nd DIR: make 2-dim leading, reshape w to p(2)-by-p(1)*m(3)
  % multiply by Q{2} and reorganize result
  w = permute(w, [2 1 3]);
  w = reshape(w,p(2),p(1)*m(3));
  w = Q{2}*w;
  w = reshape(w,m(2),p(1),m(3));
  w = permute(w,[2 1 3]);

  % 1st DIR: make 1-dim leading, reshape w to p(1)-by-m(2)*m(3)
  % multiply by Q{1} 
  w = reshape(w,p(1),m(2)*m(3));
  w = Q{1}*w;
end;

w = w(:);

function Qw = tensorProd(Q,w);
for i=1:length(Q),
  p(i,:) = size(Q{i});
end;
Qw = zeros(p(:,1)');

for i1=1:size(Qw,1),
  for i2=1:size(Qw,2),
    for i3=1:size(Qw,3),
      for j1=1:size(w,1),
        for j2=1:size(w,2),
          for j3=1:size(w,3),
            Qw(i1,i2,i3) = Q{1}(i1,j1)*Q{2}(i2,j2)*Q{3}(i3,j3)*w(j1,j2,j3);
          end;
        end;
      end;
    end;
  end;
end;

%{ 
	=======================================================================================
	FAIR: Flexible Algorithms for Image Registration, Version 2011
	Copyright (c): Jan Modersitzki
	Maria-Goeppert-Str. 1a, D-23562 Luebeck, Germany
	Email: jan.modersitzki@mic.uni-luebeck.de
	URL:   http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
	=======================================================================================
	No part of this code may be reproduced, stored in a retrieval system,
	translated, transcribed, transmitted, or distributed in any form
	or by any means, means, manual, electric, electronic, electro-magnetic,
	mechanical, chemical, optical, photocopying, recording, or otherwise,
	without the prior explicit written permission of the authors or their
	designated proxies. In no event shall the above copyright notice be
	removed or altered in any way.

	This code is provided "as is", without any warranty of any kind, either
	expressed or implied, including but not limited to, any implied warranty
	of merchantibility or fitness for any purpose. In no event will any party
	who distributed the code be liable for damages or for any claim(s) by
	any other party, including but not limited to, any lost profits, lost
	monies, lost data or data rendered inaccurate, losses sustained by
	third parties, or any other special, incidental or consequential damages
	arrising out of the use or inability to use the program, even if the
	possibility of such damages has been advised against. The entire risk
	as to the quality, the performace, and the fitness of the program for any
	particular purpose lies with the party using the code.
	=======================================================================================
	Any use of this code constitutes acceptance of the terms of the above statements
	=======================================================================================
%}