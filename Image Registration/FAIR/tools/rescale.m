function R = rescale(R,newMinR,newMaxR,gamma)
minR = min(R(:))
maxR = max(R(:))

if ~exist('newMinR'), newMinR = minR;  end;
if ~exist('newMaxR'), newMaxR = maxR;  end;
if ~exist('gamma'),   gamma   =   1;   end;

dI   = (maxR-minR); dI = dI + (dI == 0);

R = (R-minR)/dI; % R in [0,1]
R = R.^gamma;
R = (newMaxR-newMinR)*R+newMinR;

