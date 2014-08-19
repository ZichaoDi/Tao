%==============================================================================
% (c) Jan Modersitzki 2011/01/02, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% Tutorial for FAIR: density estimation using a histogram
%
%==============================================================================

clear, close all, help(mfilename);

fig = FAIRfigure(1,'position',800,'color','w');
Name    = @(n,g) sprintf('histogram-n=%s-g=%s',int2str(n),int2str(g));
FAIRprint('reset','folder',fullfile(FAIRpath,'temp'),...
  'pause',1,'obj','gcf','format','jpg','draft','off');
FAIRprint('disp');

% produce plots of different histograms
numBins = [5,10,100];
n = [10,100,1000];
for k=1:length(numBins),
  for j=1:length(n),
    T = 0.5*(sin(linspace(0,2*pi,n(j)))+1); % discretized function
    minT = 0; maxT = 1;                     % bounds for the bins
    numberBins = numBins(k);                % number of bins
    binWidth   = (maxT-minT)/numberBins;    % bin width
    bins       = 0:binWidth:maxT;           % the bins
    binsExt    = [-inf,bins(2:end-1),inf];  % don't miss anything
    rhoHat = histc(T,binsExt);
    bh = bar(bins+binWidth/2,rhoHat,0.99,'edgecolor','w','facecolor',0.8*[1,1,1]);
    axis([minT,maxT,0,inf]);
    set(gca,'fontsize',30);
    FAIRprint(Name(n(j),numberBins));
  end;
end;
