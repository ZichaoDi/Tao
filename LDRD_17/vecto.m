function vecto
TestVec
TestFor
 
function TestVec
y=zeros(1,1001);
tic;
% Vectorized form
t = 0:.01:10;
y = sin(t);
toc;
 
function TestFor
y=zeros(1,1001);
tic;
% For loop form
i = 0;
for t = 0:.01:10
i = i + 1;
y(i) = sin(t);
end
toc;

